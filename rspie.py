#!/usr/bin/env python3

import numpy as np
from subprocess import check_output
import os
import http.client, urllib
from mysecrets import *
from emfields import *
from scipy.interpolate import interp1d
from tqdm import tqdm

module_file_path = os.path.abspath(__file__)
module_directory = os.path.dirname(module_file_path)

nanolib_rectangular_template = '''

Angle = 0
Aspect = 1
Fill = 0.25
Har = 8
HarShow = 1
Offset = 0.05
Period = 0.4
PillarHeight = PillarWidth*Aspect
PillarIndex = 1.5
PillarLength = 0.6
PillarMin = 1e-06
PillarWidth = Period*Fill
SubHeight = 0.2
SubstrateHeight = 0.2
SubstrateIndex = 1.2
alpha = 0
background_index = 1
boundary_max = Period/2
boundary_max_y = Period/2
boundary_min = -Period/2
boundary_min_y = -Period/2
cad_aspectratio_x = 1
cad_aspectratio_y = 1
delta = index-background_index
dimension = 3
domain_max = PillarLength+Offset
domain_min = -Offset
eim = 0
free_space_wavelength = 0.46
height = width
index = background_index
k0 = (2*pi)/free_space_wavelength
lambda = free_space_wavelength
launch_angle = 0
launch_normalization = 2
launch_theta = 0
most_measurement_warning = 0
plot_aspectratio = 1
polarization = 0
rcwa_allhomo_warning = 0
rcwa_float = 0
rcwa_harmonics_x = Har
rcwa_harmonics_y = Har
rcwa_launch_delta_phase = 0
rcwa_launch_pol = 0
rcwa_output_absorption = 1
rcwa_output_diff_refl = 1
rcwa_output_diff_trans = 1
rcwa_output_nhx = HarShow
rcwa_output_nhy = HarShow
rcwa_output_option = 1
rcwa_output_total_refl = 1
rcwa_output_total_trans = 1
sim_tool = ST_DIFFRACTMOD
structure = STRUCT_CHANNEL
width = 1



segment 1
	structure = STRUCT_CHANNEL
	comp_name = Substrate
	begin.x = 0
	begin.z = -SubstrateHeight
	begin.height = Period
	begin.width = Period
	begin.delta = SubstrateIndex-background_index
	end.x = 0 rel begin segment 1
	end.y = 0 rel begin segment 1
	end.z = SubHeight rel begin segment 1
	end.height = Period
	end.width = Period
	end.delta = SubstrateIndex-background_index
end segment

segment 2
	comp_name = Pillar
	extended = 1
	width_taper = TAPER_LINEAR
	height_taper = TAPER_LINEAR
	begin.x = 0 rel end segment 1
	begin.y = 0 rel end segment 1
	begin.z = 0 rel end segment 1
	begin.height = max(PillarHeight,PillarMin)
	begin.width = max(PillarWidth,PillarMin)
	begin.delta = PillarIndex-background_index
	begin.euler_psi = -Angle
	end.x = 0 rel begin segment 2
	end.y = 0 rel begin segment 2
	end.z = PillarLength rel begin segment 2
	end.height = max(PillarHeight,PillarMin)
	end.width = max(PillarWidth,PillarMin)
	end.delta = PillarIndex-background_index
end segment

time_monitor 3
	profile_type = PROF_INACTIVE
	color = 2
	comp_name = Port1
	type = TIMEMON_EXTENDED
	timeaverage = 2
	monitoroutputmask = 1024
	portnum = 1
	phi = default
	begin.x = 0 rel end segment 2
	begin.y = 0 rel end segment 2
	begin.z = 0 rel end segment 2
	begin.height = Period
	begin.width = Period
end time_monitor





text_block 1
	name = MOST
	text =
RSScanOptFormat1

[MODE] 
SCAN

PREFIX mosttmp
PREFIX_STYLE 0
CLUSTER 0 0 0 0 1 ""
USERSIM_CALLSTYLE 0 0

[SIMULATION]
SIMTOOL ST_DEFAULT 
WINDOW_SIZE 0
VERBOSITY 0
PRE_WHOLE_CMD 
POST_WHOLE_CMD 
PRE_CMD 
POST_CMD 
PREPOST_ACTIVE 0
PREPOST_ERRCODES 0
EXTRA_DATAINDEX_CMDS 

[ALGORITHM]
NAME root_1d_brent
MAXSTEPS DEFAULT  1000
CONVERGENCE DEFAULT  1.0e-7

[INDEPENDENT_VARIABLES_SCAN]
IV_Declarations
SYMTAB_SCALAR Har N :  IV_LINEAR_INCR : 0 : 16 : 2 : 9 :  :  :
SYMTAB_SCALAR PillarWidth Y :  IV_LINEAR_INCR : 0 : Period : 0.02 : 21 :  :  :

[INDEPENDENT_VARIABLES_OPT]
IV_Declarations

IV_InitialValues

[MEASUREMENTS:ST_FULLWAVE]
STANDARD fw_mon_1_power_last Y 

[MEASUREMENTS:ST_DIFFRACTMOD]
STANDARD dm_de_a_total_single Y 
STANDARD dm_de_r_0_0_single Y 
STANDARD dm_de_r_total_single Y 
STANDARD dm_de_t_0_0_single Y 
STANDARD dm_de_t_total_single Y 

[MEASUREMENTS:ST_FWMPI]
STANDARD fw_mon_1_power_last Y 

[METRICS]

	end text
end text_block

text_block 2
	name = MOST_BSDFGEN
	text =
RSScanOptFormat1

[MODE] 
SCAN

PREFIX mosttmp
PREFIX_STYLE 0
CLUSTER 0 0 0 0 1 ""
USERSIM_CALLSTYLE 0 0

[SIMULATION]
SIMTOOL ST_USER bsdfgen
WINDOW_SIZE 1
VERBOSITY 0
PRE_WHOLE_CMD 
POST_WHOLE_CMD 
PRE_CMD 
POST_CMD 
PREPOST_ACTIVE 0
PREPOST_ERRCODES 0
EXTRA_DATAINDEX_CMDS 

[ALGORITHM]
NAME root_1d_brent
MAXSTEPS DEFAULT  1000
CONVERGENCE DEFAULT  1.0e-7

[INDEPENDENT_VARIABLES_SCAN]
IV_Declarations
SYMTAB_SCALAR PillarWidth Y :  IV_LINEAR_STEPS : 0 : Period : 0.02 : 21 :  :  :
SYMTAB_SCALAR PillarHeight N :  IV_LINEAR_STEPS : 0 : Period : 0.02 : 21 :  :  :
SYMTAB_SCALAR Fill N :  IV_LINEAR_STEPS : 0 : 1 : 0.05 : 21 :  :  :
SYMTAB_SCALAR Aspect N :  IV_LINEAR_STEPS : 0 : 1 : 0.25 : 5 :  :  :
SYMTAB_SCALAR Angle N :  IV_LINEAR_STEPS : 0 : 90 : 15 : 7 :  :  :

[INDEPENDENT_VARIABLES_OPT]
IV_Declarations

IV_InitialValues

[METRICS]

	end text
end text_block

'''

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
        current_dir = os.getcwd()
        self.sim_dir = os.path.join(current_dir, self.vars['filename'].replace('.ind',''))
        self.full_filename = os.path.join(self.sim_dir, self.vars['filename'])
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
        os.chdir(module_directory)
        if not os.path.exists(self.sim_dir):
            os.mkdir(self.sim_dir)
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
        if not os.path.exists(self.sim_dir):
            os.mkdir(self.sim_dir)
        os.chdir(self.sim_dir)
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
    (tuple): (x_coords, y_coords, num_array, file_format)
        file_format (str): the format of the data in the file
        x_coords (np.ndarray): 1D array of the x coordinates
        y_coords (np.ndarray): 1D array of the y coordinates
        num_array (np.ndarray): 2D array of the data with each row corresponding to strip of 
                                constant y and each column corresponding to a strip of constant x.
                                The first row corresponds to the lowest (x,y) value pair.
        
    '''
    data_lines = open(fname,'r').readlines()
    field_format = data_lines[2].split(' ')[4]
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
    if field_format == 'OUTPUT_REAL_IMAG_3D':
        num_array = num_array[:,0::2] + 1j*num_array[:,1::2]
    elif field_format == 'OUTPUT_AMP_PHASE_3D':
        phase_const = 1j/180*np.pi
        num_array = num_array[:,0::2]*np.exp(phase_const*num_array[:,1::2])
    num_array = num_array.T
    return field_format, x_coords, y_coords, num_array


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
    dtype = str(data_array.dtype)
    if 'complex' in dtype:
        fmt = '%1.5E  %1.5E  '* int(data_array.shape[0])
        format = 'OUTPUT_REAL_IMAG_3D'
    else:
        fmt = '%1.5E'
        format = 'OUTPUT_AMPLITUDE_3D'
    header = '''/rn,a,b/nx0/ls1
    /r,qa,qb
    {num_elements_x} {xmin} {xmax} 0 {format} Wavelength={wavelength}
    {num_elements_y} {ymin} {ymax}'''.format(format=format,
                                                num_elements_x=num_elements_x,
                                                num_elements_y=num_elements_y,
                                                xmin=xmin, xmax=xmax, 
                                                ymin=ymin, ymax=ymax,
                                                wavelength=wavelength)

    np.savetxt(fname,
        data_array.T,
        fmt=fmt,
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

def intervalspace(lmin, lmax, dx, symm='even'):
    '''
This  function  returns  a  numpy  array  of evenly spaced points
between  lmin  and  lmax with a spacing of dx. If symm is 'even',
the  array  will  have  an  even  number  of points symmetrically
distributed  about  the midpoint without including it. If symm is
'odd', the array will have an odd number of points with the given
points  symmetrically  distributed  about the midpoint, including
the  midpoint  itself. Given this the given array may not contain
the endpoints lmin and lmax.

    Parameters
    ----------
    lmin (float) : the minimum value of the interval
    lmax (float) : the maximum value of the interval
    dx (float)   : the spacing between points
    symm (str)   : wheter to include the midpoint or not.

    Returns
    -------
    interspace (numpy array): The evenly spaced array of points
    '''
    l = (lmax - lmin)/2
    rsteps = int(l/dx)
    if symm == 'even':
        rspace = np.linspace(dx/2, dx/2 + dx*rsteps, rsteps+1)[:-1]
        lspace = -rspace[-1::-1]
        interspace = np.concatenate((lspace, rspace))
    else:
        rspace = np.linspace(0, dx*rsteps, rsteps+1)
        lspace = -rspace[-1:0:-1]
        interspace = np.concatenate((lspace, rspace))
    interspace =  (lmin+lmax)/2 + interspace
    return interspace

def metamaker(metal_config):
    '''
    This function takes a dictionary of parameters defining a metalens
    and returns a circuit object that corresponds to it.
    The circuit also includes a launch field consisting of a plane wave
    incident at normal incidence.
    The circuit also includes a DFT monitor that monitors the field at
    the end of the pillars and saves Ex, Ey, Ez, Hx, Hy, Hz when the circuit
    is run.

    Parameters
    ----------
    metal_config : dict
        A dictionary of parameters defining a metasurface. The keys are:
        ApertureRadius : float
            The radius of the aperture in the metasurface.
        period : float
            The period of the square grid of pillars in the metasurface.
        PillarHeight : float
            The height of the pillars in the metasurface.
        circuitFname : str
            The name of the file to save the circuit to.
        PillarIndex : float
            The index of refraction of the pillars.
        SubstrateIndex : float
            The index of refraction of the substrate.
        background_index : float
            The index of refraction of the background.
        free_space_wavelength : float
            The wavelength of the light in free space.
        phase_func : function
            A function that takes in x and y coordinates and returns the
            radius that a post at that location should have. For example
            this could be the classical Fresnel phase profile.
        pillar_func : function
            A function that takes in a phase and returns the width that
            a post at that location should have to impart that phase. This
            would usually be a function that is an interpolation radii and
            phases that a previous simulation should have provided.
    Returns
    -------
    circuit : rspie.circuit.Circuit
    '''
    apertureRadius = metal_config['ApertureRadius']
    period = metal_config['period']
    phase_func = metal_config['phase_func']
    pillar_func = metal_config['pillar_func']

    # Create the pillar grid
    x = intervalspace(-apertureRadius, apertureRadius, period, 'odd')
    y = intervalspace(-apertureRadius, apertureRadius, period, 'odd')
    xgrid, ygrid = np.meshgrid(x,y)
    # calculate the required phases across the surface of the metalens
    phase_map = phase_func(xgrid, ygrid)
    # using those phases determine the radii using pillar_func
    widths = pillar_func(phase_map)

    substrate_segment = '''
    structure = STRUCT_CHANNEL
    extended = 1
    begin.x = 0
    begin.z = -(2*Offset)
    begin.height = 2*ApertureRadius
    begin.width = 2*ApertureRadius
    begin.delta = 2*SubstrateIndex - background_index
    end.x = 0 rel begin segment 1
    end.y = 0 rel begin segment 1
    end.z = 0
    end.height = 2*ApertureRadius
    end.width = 2*ApertureRadius
    end.delta = 2*SubstrateIndex-background_index
    '''
    pillar_template = '''
    extended = 1
    structure = STRUCT_CHANNEL
    position_taper = TAPER_LINEAR
    begin.x = {pillar_x}
    begin.y = {pillar_y}
    begin.z = 0
    begin.width = {pillar_diameter}
    begin.height = {pillar_diameter}
    begin.delta = PillarIndex-background_index
    end.x = {pillar_x}
    end.y = {pillar_y}
    end.z = PillarHeight
    end.width = {pillar_diameter}
    end.height = {pillar_diameter}
    end.delta = PillarIndex-background_index
    '''

    monitor_text = '''profile_type = PROF_INACTIVE
	color = 2
	type = TIMEMON_EXTENDED
	timeaverage = 2
	complexmonitor = 1
	monitoroutputmask = 0
	monitoroutputformat = OUTPUT_AMP_PHASE
	fieldoutputmask = 126
    fieldoutputformat = OUTPUT_AMP_PHASE
	frequencyanalysis = TIMEMON_FA_DFT
	dx = {output_grid_pitch}
	dy = {output_grid_pitch}
	dz = 0
	begin.x = 0
	begin.z = PillarHeight
	begin.height = ApertureRadius*2
	begin.width = ApertureRadius*2
    '''.format(**metal_config)

    config = {'vars'   : {}, 
            'segments' : [],
            'monitors' : [],
            'launch_fields' : {}
            }

    config_text = '''
    filename = {circuitFname}
    ApertureRadius = {ApertureRadius}
    Offset = lambda/2
    PillarHeight = {PillarHeight}
    PillarIndex = {PillarIndex}
    SubstrateIndex = {SubstrateIndex}
    Xmax = {ApertureRadius}
    Xmin = -{ApertureRadius}
    Ymax = {ApertureRadius}
    Ymin = -{ApertureRadius}
    Zmax = {PillarHeight} + Offset
    Zmin = -(2*Offset)
    alpha = 0
    output_grid_pitch = free_space_wavelength/5.
    background_index = {background_index}
    boundary_max = Xmax
    boundary_max_y = Ymax
    boundary_min = Xmin
    boundary_min_y = Ymin
    cad_aspectratio = 1
    delta = index-background_index
    dimension = 3
    domain_max = Zmax
    domain_min = Zmin
    eim = 0
    fdtd_display_res_auto = DISPLAY_RES_AUTO
    fdtd_monitor_time = lambda/4
    fdtd_monitor_time_auto = MONITOR_TIME_AUTO
    fdtd_pml_cells_enable = 1
    fdtd_stop_auto = 1
    fdtd_stop_time = 87
    fdtd_stop_time_auto = 1
    fdtd_time_step = 0.005681818182
    fdtd_time_step_auto = 1
    fdtd_update_time = 9*lambda/4
    fdtd_update_time_auto = DISPLAY_TIME_AUTO
    grid_size = {grid_size}
    grid_size_y = {grid_size}
    step_size = {grid_size}
    free_space_wavelength = {free_space_wavelength}
    width = 1
    height = width
    index = 1
    k0 = (2*pi)/free_space_wavelength
    lambda = free_space_wavelength
    launch_align_file = 1
    launch_height = inf
    launch_tilt = 1
    launch_type = LAUNCH_PLANEWAVE
    launch_width = inf
    sim_tool = ST_FULLWAVE
    structure = STRUCT_FIBER
    '''.format(**metal_config)

    config['vars'] = selfref_def_parser(config_text)

    pillars = []
    for (x, y, pill_width) in zip(np.ndarray.flatten(xgrid), 
                     np.ndarray.flatten(ygrid), 
                     np.ndarray.flatten(widths)):
        if x**2 + y**2 > metal_config['ApertureRadius']**2:
            continue
        pillar_x = x
        pillar_y = y
        pillars.append([pillar_x, pillar_y, pill_width])
    
    config['segments'] = [
        substrate_segment
        ]
    for idx, pillar in enumerate(pillars):
        pillar_x, pillar_y, pillar_diameter = pillar
        pillar_text = pillar_template.format(
            pillar_x = pillar_x, pillar_y = pillar_y, pillar_diameter = pillar_diameter)
        config['segments'].append(pillar_text)
    for idx, segment in enumerate(config['segments']):
        config['segments'][idx] = selfref_def_parser(segment,  config['vars'])

    config['monitors'] = [selfref_def_parser(monitor_text, config['vars'])]

    config['launch_fields'] = [selfref_def_parser('''
        launch_pathway = 0
        launch_type = LAUNCH_PLANEWAVE
        launch_tilt = 1
        launch_align_file = 1
        ''', config['vars'])
        ]
    metal = PhotoCircuit(config)
    return metal

def fresnel_profile(focal_length, medium_wavelength):
    '''
    This function returns a function that takes x and y coordinates 
    and returns the required phase pickup for the given focal length
    and medium wavelength.
    '''
    def phase_func(x, y):
        f = focal_length
        return np.mod(2*np.pi/medium_wavelength * (np.sqrt(x**2 + y**2 + f**2) - f), 2*np.pi)
    return phase_func

def atom_maker(config):
    '''
    This function determines the phase profile of meta-atoms with the given
    characteristics.
    The meta-atoms have a rectangular cross section.

    Parameters
    ----------
    config : dict with keys:
        'focal_length' : float
            focal length of the lens
        'medium_wavelength' : float
            wavelength of the medium
        'Period' : float
            period of the square grid used to simulate the phase response of meta-atoms
        'free_space_wavelength' : float
            wavelength of the free space
        'PillarIndex' : float
            index of the pillar
        'SubstrateIndex' : float
            index of the substrate
        'background_index' : float
            index of the background
        'Angle' : float
            angle of the rectangular meta-atom
        'fill_steps' : int
            number of fill_factors to use
        'Har' : int
            how many harmonics are used in the RCWA simulation
        'Aspect' : float
            aspect ratio of the rectangular meta-atom
        'min_width' : float
            minimum width of the pillars
        'max_width' : float
            maximum width of the pillars
    
    Returns
    -------
    overlap_phases : np.array
        overlap phases for each fill factor, given in radians
    pillar_widths : list of floats
        widths of the pillars corresponding to the phases
    '''
    os.chdir(module_directory)
    
    overlap_magnitudes = []
    overlap_phases     = []

    fill_steps = config['fill_steps']
    min_fill = config['min_width']/config['Period']
    max_fill = config['max_width']/config['Period']
    fills = np.linspace(min_fill, max_fill, fill_steps)

    pillar_widths = fills * config['Period']

    if not os.path.exists('nanolib_rectangular_template.ind'):
        print("nanolib_rectangular_template.ind not found. Creating...")
        open('nanolib_rectangular_template.ind', 'w').write(nanolib_rectangular_template)

    for fill in tqdm(fills):
        simulscript = 'nanolib_rectangular_template.ind'
        config['Fill'] = fill
        omag, ophase = rectangular_meta_atom(config, 
                                             simulscript, 
                                             hide=True, 
                                             cleanup=True)
        overlap_magnitudes.append(omag)
        overlap_phases.append(ophase)
    overlap_magnitudes = np.array(overlap_magnitudes)
    overlap_phases = np.array(overlap_phases)
    try:
        send_message('finished!')
    except:
        print("Error sending message.")
    overlap_phases = overlap_phases/360*(2*np.pi)
    overlap_phases = np.unwrap(overlap_phases)
    overlap_phases = overlap_phases - np.min(overlap_phases)
    return overlap_phases, pillar_widths
    
