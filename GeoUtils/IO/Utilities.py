__doc__ = 'Utility I/O functions'

import json

def topologystring(geotype, name, patches, entries):
  xml = '  <set name="%s" type="%s">\n' % (name, geotype)
  for patch in patches:
    xml += '    <item patch="%i">%s</item>\n' % (patch, entries)
  xml += '  </set>\n'
  return xml

def topologystringfrompairs(geotype, name, entries):
  xml = '  <set name="%s" type="%s">\n' % (name, geotype)
  for (patch, face) in entries:
    xml += '    <item patch="%i">%s</item>\n' % (patch, face)
  xml += '  </set>\n'
  return xml

def writeAsc(X, U, fname):
  f = open(fname,'w')
  for i in range(0,len(U)):
    f.write('%f %f\n' % (X[i], U[i]))
  f.close()

def ParseArgs(args, defaults):
  """ Parses command line arguments in a standardized way. Takes a list of command line
      arguments and a dictionary mapping arguments to default values. The defaults dict
      will be modified. The arguments must be of the form: 'argname=newvalue'

      This works for string, integer and float values (as determined by the default).
      For boolean values, an argument 'argname' will set argname to True, while an argument
      'noargname' will set it to False.

      Arguments are parsed in order, so later values can override earlier ones.

      A special argument of form 'paramfile=somefile.json' will load 'somefile.json' and apply
      the arguments in that file, in order. To specify boolean arguments in JSON, give true
      or false as the value.

      @param args: Command-line arguments
      @type args: List of String
      @param defaults: Defaults
      @type defaults: Dict
      @returns: None
      @rtype: NoneType
  """

  def set_def(key, value):
    for k, v in defaults.iteritems():
      if value is None and key in [k, 'no' + k] and isinstance(v, bool):
        defaults[k] = k == key
      elif k == key:
        if isinstance(v, bool):
          value = {'True':True, 'False':False}.get( value, value ) # try string conversion
          defaults[k] = bool(value)
        elif isinstance(v, float):
          defaults[k] = float(value)
        elif isinstance(v, int):
          defaults[k] = int(value)
        else:
          defaults[k] = str(value)

  for arg in args:
    if arg.startswith('paramfile=') :
      with open(arg[10:], 'r') as f:
        dc = json.load(f)
      for k, v in dc.iteritems():
        set_def(k, v)

    split = arg.split('=', 1)
    set_def(split[0], split[1] if len(split) > 1 else None)

def readMCfile(filename,remlast=True, silent=False):
    """
    Read text file with multiple columns
    If remlast = True, the last row will be deleted
    Return as tuple (x,y)
    """

    # Open file for reading
    try:
        readfile = open(filename,'r')
    except IOError:
        if silent != True:
            print 'Utils.FileInputOutput.py::readMCfile: Could not open file ' + filename + ' for reading.'
        #sys.exit(1)
        return (None, None)

    try:
        data = scitools.filetable.read(readfile)
    except Exception, e:
        print 'Utils.FileInputOutput.py::readMCfile: Could not read from file ' + filename + '.'
        print str(e)
        data = None

        #sys.exit(1)

    readfile.close()
    if data == None:
        print 'Utils.FileInputOutput.py::readMCfile: Trying alterative file reading method for' + filename + '.'
        data =  readLinesWithSeparation(filename,remlast=True, silent=silent)
        #pdb.set_trace()
    try:
      if remlast == True:
          x = data[0:-1,0]        # Remove last line which often is problematic
          y = data[0:-1,1:]       # Remove last line which often is problematic
      else:
          x = data[:,0]
          y = data[:,1:]
      x.reshape(x.shape[0],1)
    except TypeError:
      return (None, None)

    return (x,y)
