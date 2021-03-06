#!/usr/bin/env python
# encoding: utf-8
'''Convert the binary output of DynEarthSol3D to VTK files.

usage: 2vtk.py [-a -c -h] modelname [start [end [delta]]]]

options:
    -a          save data in ASCII format (default: binary)
    -c          save files in current directory (default: same directory as
                the data files)
    -h,--help   show this help
'''

from __future__ import print_function, unicode_literals
import sys, os
import base64, zlib
import numpy as np
from Dynearthsol import Dynearthsol

# Save in ASCII or encoded binary.
# Some old VTK programs cannot read binary VTK files.
output_in_binary = True

# Save the resultant vtu files in current directory?
output_in_cwd = False

# Save indivisual components?
output_tensor_components = False

# Save P-/T-axis
output_pt_axis = False

# Save markers?
output_markers = False


def main(modelname, start, end, delta):
    prefix = modelname
    if output_in_cwd:
        output_prefix = os.path.basename(modelname)
    else:
        output_prefix = modelname

    des = Dynearthsol(modelname)

    if end == -1:
        end = len(des.frames)

    for i in range(start, end, delta):
        frame = des.frames[i]
        nnode = des.nnode_list[i]
        nelem = des.nelem_list[i]
        step = des.steps[i]
        time_in_yr = des.time[i] / (365.2425 * 86400)

        des.read_header(frame)
        suffix = '{0:0=6}'.format(frame)
        print('Converting frame #{0}'.format(suffix), end='\r', file=sys.stderr)

        filename = '{0}.{1}.vtu'.format(output_prefix, suffix)
        fvtu = open(filename, 'wb')

        try:
            vtu_header(fvtu, nnode, nelem, time_in_yr, step)

            #
            # node-based field
            #
            fvtu.write(b'  <PointData>\n')

            convert_field(des, frame, 'temperature', fvtu)
            #convert_field(des, frame, 'z0', fvtu)
            #convert_field(des, frame, 'bcflag', fvtu)
            #convert_field(des, frame, 'mass', fvtu)
            #convert_field(des, frame, 'tmass', fvtu)
            #convert_field(des, frame, 'volume_n', fvtu)

            # node number for debugging
            vtk_dataarray(fvtu, np.arange(nnode, dtype=np.int32), 'node number')

            fvtu.write(b'  </PointData>\n')

            #
            # node coordinate
            #
            fvtu.write(b'  <Points>\n')
            convert_field(des, frame, 'coordinate', fvtu)
            fvtu.write(b'  </Points>\n')

            #
            # element connectivity & types
            #
            fvtu.write(b'  <Cells>\n')
            convert_field(des, frame, 'connectivity', fvtu)
            vtk_dataarray(fvtu, (des.ndims+1)*np.array(range(1, nelem+1), dtype=np.int32), 'offsets')
            if des.ndims == 2:
                # VTK_ TRIANGLE == 5
                celltype = 5
            else:
                # VTK_ TETRA == 10
                celltype = 10
            vtk_dataarray(fvtu, celltype*np.ones((nelem,), dtype=np.int32), 'types')
            fvtu.write(b'  </Cells>\n')

            vtu_footer(fvtu)
            fvtu.close()

        except:
            # delete partial vtu file
            fvtu.close()
            os.remove(filename)
            raise

    print()
    return


def output_vtp_file(des, frame, filename, markersetname, time_in_yr, step):
    fvtp = open(filename, 'wb')

    class MarkerSizeError(RuntimeError):
        pass

    try:
        # read data
        marker_data = des.read_markers(frame, markersetname)
        nmarkers = marker_data['size']

        if nmarkers <= 0:
            raise MarkerSizeError()

        # write vtp header
        vtp_header(fvtp, nmarkers, time_in_yr, step)

        # point-based data
        fvtp.write('  <PointData>\n')
        name = markersetname + '.mattype'
        vtk_dataarray(fvtp, marker_data[name], name)
        name = markersetname + '.elem'
        vtk_dataarray(fvtp, marker_data[name], name)
        name = markersetname + '.id'
        vtk_dataarray(fvtp, marker_data[name], name)
        fvtp.write('  </PointData>\n')

        # point coordinates
        fvtp.write('  <Points>\n')
        field = marker_data[markersetname + '.coord']
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            tmp = np.zeros((nmarkers, 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        vtk_dataarray(fvtp, tmp, markersetname + '.coord', 3)
        fvtp.write('  </Points>\n')

        vtp_footer(fvtp)
        fvtp.close()

    except MarkerSizeError:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        # skip this frame

    except:
        # delete partial vtp file
        fvtp.close()
        os.remove(filename)
        raise

    return


def convert_field(des, frame, name, fvtu):
    field = des.read_field(frame, name)
    if name in ('coordinate', 'velocity', 'velocity averaged', 'force'):
        if des.ndims == 2:
            # VTK requires vector field (velocity, coordinate) has 3 components.
            # Allocating a 3-vector tmp array for VTK data output.
            i = des.frames.index(frame)
            tmp = np.zeros((des.nnode_list[i], 3), dtype=field.dtype)
            tmp[:,:des.ndims] = field
        else:
            tmp = field

        # Rename 'velocity averaged' to 'velocity'
        if name == 'velocity averaged': name = 'velocity'

        vtk_dataarray(fvtu, tmp, name, 3)
    else:
        vtk_dataarray(fvtu, field, name)
    return


def vtk_dataarray(f, data, data_name=None, data_comps=None):
    if data.dtype in (np.int32,):
        dtype = 'Int32'
    elif data.dtype in (np.single, np.float32):
        dtype = 'Float32'
    elif data.dtype in (np.double, np.float64):
        dtype = 'Float64'
    else:
        raise Error('Unknown data type: ' + name)

    name = ''
    if data_name:
        name = 'Name="{0}"'.format(data_name)

    ncomp = ''
    if data_comps:
        ncomp = 'NumberOfComponents="{0}"'.format(data_comps)

    if output_in_binary:
        fmt = 'binary'
    else:
        fmt = 'ascii'
    header = '<DataArray type="{0}" {1} {2} format="{3}">\n'.format(
        dtype, name, ncomp, fmt)
    f.write(header.encode('ascii'))
    if output_in_binary:
        header = np.zeros(4, dtype=np.int32)
        header[0] = 1
        a = data.tostring()
        header[1] = len(a)
        header[2] = len(a)
        b = zlib.compress(a)
        header[3] = len(b)
        f.write(base64.standard_b64encode(header.tostring()))
        f.write(base64.standard_b64encode(b))
    else:
        data.tofile(f, sep=b' ')
    f.write(b'\n</DataArray>\n')
    return


def vtu_header(f, nnode, nelem, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<UnstructuredGrid>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {3}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}" NumberOfCells="{1}">
'''.format(nnode, nelem, time, step).encode('ascii'))
    return


def vtu_footer(f):
    f.write(
b'''</Piece>
</UnstructuredGrid>
</VTKFile>
''')
    return


def vtp_header(f, nmarkers, time, step):
    f.write(
'''<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">
<PolyData>
<FieldData>
  <DataArray type="Float32" Name="TIME" NumberOfTuples="1" format="ascii">
    {1}
  </DataArray>
  <DataArray type="Float32" Name="CYCLE" NumberOfTuples="1" format="ascii">
    {2}
  </DataArray>
</FieldData>
<Piece NumberOfPoints="{0}">
'''.format(nmarkers, time, step).encode('ascii'))
    return


def vtp_footer(f):
    f.write(
b'''</Piece>
</PolyData>
</VTKFile>
''')
    return


def first_invariant(t):
    nstr = t.shape[1]
    ndims = 2 if (nstr == 3) else 3
    return np.sum(t[:,:ndims], axis=1) / ndims


def second_invariant(t):
    '''The second invariant of the deviatoric part of a symmetric tensor t,
    where t[:,0:ndims] are the diagonal components;
      and t[:,ndims:] are the off-diagonal components.'''
    nstr = t.shape[1]

    # second invariant: sqrt(0.5 * t_ij**2)
    if nstr == 3:  # 2D
        return np.sqrt(0.25 * (t[:,0] - t[:,1])**2 + t[:,2]**2)
    else:  # 3D
        a = (t[:,0] + t[:,1] + t[:,2]) / 3
        return np.sqrt( 0.5 * ((t[:,0] - a)**2 + (t[:,1] - a)**2 + (t[:,2] - a)**2) +
                        t[:,3]**2 + t[:,4]**2 + t[:,5]**2)


def compute_pt_axis(stress):
    '''The compression and dilation axis of the deviatoric stress tensor.'''

    nelem = stress.shape[0]
    nstr = stress.shape[1]
    # VTK requires vector field (velocity, coordinate) has 3 components.
    # Allocating a 3-vector tmp array for VTK data output.
    p_axis = np.zeros((nelem, 3), dtype=stress.dtype)
    t_axis = np.zeros((nelem, 3), dtype=stress.dtype)

    if nstr == 3:  # 2D
        sxx, szz, sxz = stress[:,0], stress[:,1], stress[:,2]
        mag = np.sqrt(0.25*(sxx - szz)**2 + sxz**2)
        theta = 0.5 * np.arctan2(2*sxz, sxx-szz)
        cost = np.cos(theta)
        sint = np.sin(theta)

        p_axis[:,0] = mag * sint
        p_axis[:,1] = mag * cost
        t_axis[:,0] = mag * cost
        t_axis[:,1] = -mag * sint

    else:  # 3D
        # lower part of symmetric stress tensor
        s = np.zeros((nelem, 3,3), dtype=stress.dtype)
        s[:,0,0] = stress[:,0]
        s[:,1,1] = stress[:,1]
        s[:,2,2] = stress[:,2]
        s[:,1,0] = stress[:,3]
        s[:,2,0] = stress[:,4]
        s[:,2,1] = stress[:,5]

        # eigenvalues and eigenvectors
        if np.version.version[:3] >= '1.8':
            # Numpy-1.8 or newer
            w, v = eigh(s)
        else:
            # Numpy-1.7 or older
            w = np.zeros((nelem,3), dtype=stress.dtype)
            v = np.zeros((nelem,3,3), dtype=stress.dtype)
            for e in range(nelem):
                w[e,:], v[e,:,:] = eigh(s[e])

        # isotropic part to be removed
        m = np.sum(w, axis=1)

        p = w.argmin(axis=1)
        t = w.argmax(axis=1)
        #print(w.shape, v.shape, p.shape)

        for e in range(nelem):
            p_axis[e,:] = (w[e,p[e]] - m[e]) * v[e,:,p[e]]
            t_axis[e,:] = (w[e,t[e]] - m[e]) * v[e,:,t[e]]

    return p_axis, t_axis


if __name__ == '__main__':

    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            if arg.lower() in ('-h', '--help'):
                print(__doc__)
                sys.exit(0)

    if '-a' in sys.argv:
        output_in_binary = False
    if '-c' in sys.argv:
        output_in_cwd = True
    if '-p' in sys.argv:
        output_pt_axis = True
    if '-t' in sys.argv:
        output_tensor_components = True
    if '-m' in sys.argv:
        output_markers = True

    # delete options
    for i in range(len(sys.argv)):
        if sys.argv[1][0] == '-':
            del sys.argv[1]
        else:
            # the rest of argv cannot be options
            break

    modelname = sys.argv[1]

    if len(sys.argv) < 3:
        start = 0
    else:
        start = int(sys.argv[2])

    if len(sys.argv) < 4 or int(sys.argv[3]) == -1:
        end = -1
    else:
        end = int(sys.argv[3]) + 1

    if len(sys.argv) < 5:
        delta = 1
    else:
        delta = int(sys.argv[4])

    main(modelname, start, end, delta)
