#!/usr/bin/env python3

import numpy as np
from subprocess import check_output
import os
import http.client, urllib
from mysecrets import *
from emfields import *

def send_message(message):
    app_token = pushover_token
    conn = http.client.HTTPSConnection("api.pushover.net",443)
    endpoint = "/1/messages.json"
    conn.request("POST", endpoint,
      urllib.parse.urlencode({
        "token": app_token,
        "user": pushover_user,
        "message": message,
      }), { "Content-type": "application/x-www-form-urlencoded" })
    return conn.getresponse().read().decode()

def rectangular_meta_atom(config, simulscript = '', hide=True, cleanup=True):
    '''
    This  function  runs  a  simulation  using DiffractMOD to estimate the
    phase  difference  that  a  rectangular  pillar  of given geometry and
    composition will provide to a plane wave at normal incidence.

    This is done by simulating a periodic grating in a square lattice with
    the pillar in the center.

    The launch field, as configured in the auxiliary script, is assumed to
    be  a  plane wave at normal incidence, and is launched from inside the
    substrate. The launch field has linear polarization along the x-axis.

    The  pillar  is  assumed  to be on top of a substrate whose refractive
    index can also be provided.

    All the refractive indices are assumed to be real and isotropic.

    The width is the length of the pillar along the x-axis, and the height
    is  the  length  along the y-axis. This before the rotation defined by
    Angle is applied.

    Parameters
    ----------
    config (dict) : with the following keys (not all which need to be provided)
        'Angle' (float)  : angle (in degrees) that orients the pillar (default=0)
        'Aspect' (float) : width/height of the pillar (default=1)
        'Fill'   (float) : width/Period (used to fix width) (default=0.25)
        'Period' (float) : period of the square lattice in um.
        'PillarHeight' (float): length of pillar (in um) in the z-dir (measured from the top of the substrate)
        'PillarIndex'  (float): refractive index of the pillar
        'SubstrateIndex' (float): refractive index of the substrate
        'background_index' (float): refractive index of the background
    simulscript (str): path to the template circuit to run the simulation
    hide (bool)      : whether to hide the simulation window or not
    cleanup (bool)   : whether to delete the output files or not

    Returns
    -------
    (float, float) : (overlap_magnitude, overlap_phase_in_degrees)

    '''

    # Put together the command to run the simulation
    if simulscript == '':
        simulscript = 'nanolib_rectangular_template.ind'
    if hide:
        cmd = ['dfmod', '-hide', simulscript, 'wait=0']
    else:
        cmd = ['dfmod', simulscript, 'wait=0']
    int_params = ['%s=%f' % (k, v) for k, v in config.items() if type(v) in [int]]
    num_params = ['%s=%f' % (k, v) for k, v in config.items() if type(v) in [float, np.float64]]
    str_params = ['%s=%s' % (k, v) for k, v in config.items() if type(v) == str]
    cmd = cmd + int_params + num_params + str_params

    # Run it
    check_output(cmd, shell=True)

    # read the results from the log file
    log = open('%s.txt' % config['prefix'],'r').read()
    line = log.split('\n')[-2]
    overlap_mag = float(line.split(',')[0].split('=')[-1].strip())
    overlap_phase_in_degrees = float(line.split(',')[-1].split('=')[1].strip())

    if cleanup:
        extra_files = [f for f in os.listdir() if f.startswith(config['prefix'])]
        for f in extra_files:
            os.remove(f)
    return (overlap_mag, overlap_phase_in_degrees)

def loadfld(filename):
    '''
    Parameters
    ----------
    filename (str): path to the .fld file

    Returns
    -------
    extent, data (tuple): (extent, data)
        extent (tuple): (xmin, xmax, zmin, zmax)
        data (np.ndarray): 2D array of the data
    '''
    with open(filename, 'r') as f:
        lines = f.readlines()
        metadata = lines[:4]
        xmin, xmax = metadata[2].split(' ')[1:3]
        zmin, zmax = metadata[3].split(' ')[1:3]
        lines = lines[4:]
        lines = [line.strip() for line in lines]
        lines = [line.split() for line in lines]
        lines = [[float(x) for x in line] for line in lines]
        extent = (float(xmin), float(xmax), float(zmin), float(zmax))
        return extent, np.array(lines).T

cmd_tool_dict = {
    'ST_BEAMPROP': 'bsimw32',
    'ST_GRATINGMOD': 'grmod',
    'ST_DIFFRACTMOD': 'dfmod',
    'ST_MODEPROP': 'bsimw32',
    'ST_FULLWAVE': 'fullwave',
    'ST_BANDSOLVE': 'bandsolve',
    'ST_FEMSIM': 'femsim'
    }

class PhotoCircuit():
    '''
    Helper class to create a circuit file for RSoft.

    To create a circuit file you need to provide a configuration
    dictionary with the following keys:

        'vars':  a dictionary with keys and values that will be used
        to   form  the  symbols  table  of  the  circuit  file.  The
        dictionary  must  include  the  key 'circuit_filename' which
        will be used to export the circuit file to disk.

        'segments':  a  list of dictionaries with keys sufficient to
        define   them,  e.g.:  'structure',  'material',  'begin.x',
        'begin.z',  'begin.width', 'begin.height', 'end.x', 'end.y',
        'end.z',  'position_taper'.  To  figure  out  which keys are
        necessary  the  best  way  is to create a similar circuit in
        RSoft  and  inspect the values that it uses to specify it in
        the corresponding circuit file.

        'monitors':  a  list of dictionaries with keys sufficient to
        define  them.  Again  best  way to figure out which keys are
        sufficient   is  to  create  a  similar  circuit  in  RSoft.
        'launch_fields': a list of dictionaries with keys sufficient
        to  define them. Again best way to figure out which keys are
        sufficient is to create a similar circuit in RSoft.
    
    This class has the following function attributes:
        parse_config:  parses  the configuration dictionary and sets
        the corresponding attributes of the class.

        block_parser:  parses  the  blocks  that define the circuit.
        make_circuit_text:  creates  the  circuit file text from the
        configuration   dictionary   and   the   attributes  set  by
        parse_config.

        parse_vars:  parses  the  variables  dictionary and sets the
        corresponding  attributes  of the class. write_circuit_file:
        writes  the circuit file to disk. run: runs the circuit file
        using  the  simulation  tool  specified in the configuration
        dictionary.
    '''
    def __init__(self, config):
        self.vars = config['vars']
        self.full_filename = os.path.join(os.getcwd(), self.vars['filename'])
        self.segments = config['segments']
        self.monitors = config['monitors']
        self.launch_fields = config['launch_fields']
        self.parse_config()
        self.circuit_text = self.make_circuit_text()
        self.executable = cmd_tool_dict[self.vars['sim_tool']]
    
    def parse_config(self):
        '''
        Put everything together
        '''
        self.parse_vars()
        self.num_segments, self.segment_block = self.block_parser(self.segments, 'segment')
        self.num_monitors, self.monitor_block = self.block_parser(self.monitors, 'time_monitor', self.num_segments)
        self.num_launch_fields, self.launch_field_block = self.block_parser(self.launch_fields, 'launch_field')

    def block_parser(self, block_dicts, block_header, offset=0):
        '''
        Create  a string of definitions using the key-value pairs in
        the  given  dictionary,  prepending  and  appending adequate
        headers and footers.

        Parameters
        ----------
        block_dicts (list): a list with dictionaries defining the block
        block_header (str): header of the block (e.g. 'segment')
        offset (int)      : offset to add to the index of the block

        Returns
        -------
        (tuple) (block_length (int) , block(str))
        '''
        blocks = []
        for idx, chunk in enumerate(block_dicts):
            block = ['%s %d' % (block_header, idx+1+offset)]
            for var, var_value in chunk.items():
                if type(var_value) == str:
                    block.append('\t%s = %s' % (var, var_value))
                elif type(var_value) == int:
                    block.append('\t%s = %d' % (var, var_value))
                else:
                    block.append('\t%s = %f' % (var, var_value))
            block.append('end %s' % block_header)
            blocks.append('\n'.join(block))
        return len(blocks), '\n\n'.join(blocks) 
    
    def make_circuit_text(self):
        '''
        Put together the circuit text from the its segments,
        monitors, launch fields and variables.

        Parameters
        ----------
        None

        Returns
        -------
        circuit_text (str): the circuit text
        '''
        return '\n\n'.join([self.var_block, self.segment_block, self.monitor_block, self.launch_field_block])

    def parse_vars(self):
        '''
        Parse self.vars as instance attributes.

        Parameters
        ----------
        None
        Returns
        -------
        None
        '''
        if type(self.vars) == 'str':
            self.var_block = self.vars
        else:
            var_block = []
            for var, var_value in self.vars.items():
                if type(var_value) == str:
                    var_block.append('%s = %s' % (var, var_value))
                elif type(var_value) == int:
                    var_block.append('%s = %d' % (var, var_value))
                else:
                    var_block.append('%s = %f' % (var, var_value))
            self.var_block = '\n'.join(var_block)
    
    def save_to_file(self):
        '''
        Save circuit to self.full_filename.

        Returns
        -------
        None
        '''
        with open(self.full_filename, 'w') as f:
            f.write(self.circuit_text)
    
    def run(self):
        '''
        Run the circuit.

        Parameters
        ----------
        None

        Returns
        -------
        None
        '''
        # save the file to disk
        self.save_to_file()
        # run the circuit
        # in the philosophy of this script, any change in parameters
        # should be reflected in the circuit file, so we don't need to
        # use the options to change parameter values.
        # If using command line flags is necessary, a custom
        # command line parser should be written.
        cmd = [self.executable, self.full_filename, 'wait=0']
        check_output(cmd, shell=True)

    def open_in_RSoft(self):
        '''
        Opens the created circuit in RSoft.
        '''
        self.save_to_file()
        os.startfile(self.full_filename)

    def __str__(self):
        return self.circuit_text
    
    def __repr__(self):
        return self.circuit_text

def selfref_def_parser(def_text, aux_vars={}, max_iteration_depth=10):
    """
    Give a text string of definitions, parse them into a dictionary.
    Each line in the text string should be of the form: var_name = var_value
    where var_name should be a valid Python variable name, and var_value
    a string or an expression that may depend on the other variables defined
    in the given text.

    The given text might have commentary lines starting with #, which will
    make them be omitted. The text may also have empty lines.

    If provided aux_vars is a dictionary whose values and keys will also be used
    to make sense of the definitions in the given text.
    
    Parameters
    ----------
    def_text (str): text string of definitions
    aux_vars (dict): auxiliary variables to be used in the expressions
    max_iteration_depth (int): number of iterations to try to parse the definitions

    Returns
    -------
    vars_dict (dict): dictionary of definitions

    Example
    -------
    >>> def_text = '''x = 1
    y = x+1
    '''
    >>> def_dict = recursive_def_parser(def_text)
    >>> print(def_dict)
        {'x': 1, 'y': 2}
    """
    vars_dict = {}
    varlines = [s for s in def_text.split('\n') if s.strip() != '']
    varlines = [s for s in varlines if s[0] != '#']
    old_fails = -1
    exit_next = False
    for iter in range(max_iteration_depth):
        fails = 0
        for line in varlines:
            var_name, var_val = line.split('=')
            var_name, var_val = var_name.strip(), var_val.strip()
            try:
                vars_dict[var_name] = eval(var_val, vars_dict, aux_vars)
            except:
                vars_dict[var_name] = var_val
                fails += 1
        del(vars_dict['__builtins__'])
        if exit_next:
            break
        if fails == old_fails:
            exit_next = True
        old_fails = fails
    return vars_dict

def load_3d_dat(fname):
    '''
    This function can be used to load a 3D .dat file from RSoft.
    It assumes that the file has real-valued data.
    Parameters
    ----------
    fname (str): path to the .dat file
    Returns
    -------
    x_coords, y_coords, z_coords, num_array (tuple): (x_coords, y_coords, z_coords, num_array)
        x_coords (np.ndarray): 1D array of the x coordinates
        y_coords (np.ndarray): 1D array of the y coordinates
        z_coords (np.ndarray): 1D array of the z coordinates
        num_array (np.ndarray): 3D array of the data
    '''
    data_lines = open(fname,'r').readlines()
    metadata_X = data_lines[2].split(' ')[:3]
    metadata_Y = data_lines[3].split(' ')[:3]
    metadata_Z = data_lines[4].split(' ')[:3]
    x_data_points, x_min, x_max = [float(x) for x in metadata_X]
    x_data_points = int(x_data_points)
    y_data_points, y_min, y_max = [float(x) for x in metadata_Y]
    y_data_points = int(y_data_points)
    z_data_points, z_min, z_max = [float(x) for x in metadata_Z]
    z_data_points = int(z_data_points)
    data_lines = data_lines[5:]
    num_array = []
    num_plane = []
    for line in data_lines:
        nums = line.split('  ')
        nums = [float(x) for x in nums]
        num_plane.append(nums)
        if len(num_plane) % x_data_points == 0:
            num_array.append(num_plane)
            num_plane = []
    num_array = np.array(num_array)
    x_coords = np.linspace(x_min, x_max, x_data_points)
    y_coords = np.linspace(y_min, y_max, y_data_points)
    z_coords = np.linspace(z_min, z_max, z_data_points)
    num_array = np.transpose(num_array, (1,2,0))
    return x_coords, y_coords, z_coords, num_array

def load_2d_dat(fname):
    '''
    This function can be used to load a 2D .dat file from RSoft.
    It assumes that the file has real-valued data.

    Parameters
    ----------
    fname (str): path to the .dat file

    Returns
    -------
    x_coords, y_coords, z_coords, num_array (tuple): (x_coords, y_coords, z_coords, num_array)
        x_coords (np.ndarray): 1D array of the x coordinates
        y_coords (np.ndarray): 1D array of the y coordinates
        num_array (np.ndarray): 2D array of the data with each row corresponding to strip of 
                                constant y and each column corresponding to a strip of constant x.
                                The first row corresponds to the lowest (x,y) value pair.
    '''
    data_lines = open(fname,'r').readlines()
    metadata_X = data_lines[2].split(' ')[:3]
    metadata_Y = data_lines[3].split(' ')[:3]
    x_data_points, x_min, x_max = [float(x) for x in metadata_X]
    x_data_points = int(x_data_points)
    y_data_points, y_min, y_max = [float(x) for x in metadata_Y]
    y_data_points = int(y_data_points)
    data_lines = data_lines[4:]
    num_array = []
    for line in data_lines:
        nums = line.split('  ')
        nums = [float(x) for x in nums]
        num_array.append(nums)
    num_array = np.array(num_array)
    x_coords = np.linspace(x_min, x_max, x_data_points)
    y_coords = np.linspace(y_min, y_max, y_data_points)
    num_array = num_array.T
    return x_coords, y_coords, num_array

def save_2D_array_to_dat(fname, data_array, wavelength, xmin, xmax, ymin, ymax):
    '''
    This function saves a 2D array to a dat file that can be read by RSoft.
    xmin, xmax, ymin, ymax need to match the bounds of the simulation volume
    in the RSoft simulation where the file will be used as input.

    Parameters
    ----------
    fname (str)  : The filename where the data will be saved.
    data_array (numpy.ndarray) : The 2D array that will be saved.
    wavelength (float) : The wavelength of the simulation.
    xmin (float) : The minimum x value of the simulation.
    xmax (float) : The maximum x value of the simulation.
    ymin (float) : The minimum y value of the simulation.
    ymax (float) : The maximum y value of the simulation.

    Returns
    -------
    None
    '''
    num_elements_y, num_elements_x = data_array.shape
    header = '''/rn,a,b/nx0/ls1
/r,qa,qb
{num_elements_x} {xmin} {xmax} 0 OUTPUT_AMPLITUDE_3D Wavelength={wavelength}
{num_elements_y} {ymin} {ymax}'''.format(num_elements_x=num_elements_x,
                                             num_elements_y=num_elements_y,
                                             xmin=xmin, xmax=xmax, 
                                             ymin=ymin, ymax=ymax,
                                             wavelength=wavelength)
    np.savetxt(fname,
        data_array.T,
        fmt='%1.5E',
        delimiter='  ',
        newline='\n',
        header=header,
        footer='',
        comments='')
    return None

# dipole fields

def dipole_field(x,y,z,xd,yd,zd,thetadip,phidip,omega):
    '''
    Calculates the electric and magnetic fields of an electric dipole.

    Parameters
    ----------
    x, y, z          : (float, float, float) coordinates of the point where the field is calculated
    xd, yd, zd       : (float, float, float) coordinates of the dipole
    thetadip, phidip : (float, float) polar and azimuthal angles of the dipole moment
    omega            : (float) frequency of the dipole

    Returns
    -------
    EBfield : (np.array) electric and magnetic fields (Ex, Ey, Ez, Bx, By, Bz)

    '''
    EBfield = np.array([Edipx(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Edipy(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Edipz(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipx(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipy(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipz(x,xd,y,yd,z,zd,thetadip,phidip,omega)])
    return EBfield

def dipole_field_far(x,y,z,xd,yd,zd,thetadip,phidip,omega):
    '''
    Calculates the electric and magnetic fields of an electric dipole in the
    radiation zone.

    Parameters
    ----------
    x, y, z          : (float, float, float) coordinates of the point where the field is calculated
    xd, yd, zd       : (float, float, float) coordinates of the dipole
    thetadip, phidip : (float, float) polar and azimuthal angles of the dipole moment
    omega            : (float) frequency of the dipole

    Returns
    -------
    EBfield : (np.array) electric and magnetic fields (Ex, Ey, Ez, Bx, By, Bz)

    '''
    EBfield = np.array([Edipfarx(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Edipfary(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Edipfarz(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipfarx(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipfary(x,xd,y,yd,z,zd,thetadip,phidip,omega),
                        Bdipfarz(x,xd,y,yd,z,zd,thetadip,phidip,omega)])
    return EBfield

# dipField = np.vectorize(dipField, excluded=['z', 'xd','yd','zd','thetadip','phidip','omega'])

