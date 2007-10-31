import glob
import os

Import("*")

# Note that this script is executed from within the <project>/build/<something>/ folder
curdir = os.getcwd()

# Get our project directory name
project_name = os.path.basename(os.path.dirname(os.path.dirname(curdir)))

#print curdir, project_name

# Determine Chaste libraries to link against.
# Note that order does matter!
chaste_libs = [project_name] + comp_deps['core']
#chaste_libs = [project_name, 'cancer'] + comp_deps['cancer']
#chaste_libs = [project_name, 'heart'] + comp_deps['heart']
#chaste_libs = [project_name, 'cancer', 'heart'] + comp_deps['heart']


# Look for .cpp files within the src folder
os.chdir('../..') # This is so .o files are built in <project>/build/<something>/
files, extra_cpppath = SConsTools.FindSourceFiles('src')
extra_cpppath.append('src')

# Look for source files that tests depend on under test/.
testsource, test_cpppath = SConsTools.FindSourceFiles('test', ['data'])
extra_cpppath.extend(test_cpppath)
del test_cpppath

os.chdir(curdir)

# Look for files containing a test suite
# A list of test suites to run will be found in a test/<name>TestPack.txt
# file, one per line.
# Alternatively, a single test suite may have been specified on the command
# line.
test_this_comp = False
for targ in BUILD_TARGETS:
    if str(targ) in ['projects/'+project_name, '.', Dir('#').abspath]:
        test_this_comp = True
testfiles = set()
if single_test_suite:
  if single_test_suite_dir == project_name:
    testfiles.add(single_test_suite)
    # Remove any old test output file to force a re-run
    try:
      os.remove(single_test_suite[:-4] + '.log')
    except OSError:
      pass
elif test_this_comp:
  packfiles = []
  if all_tests:
    for packfile in glob.glob('../../test/*TestPack.txt'):
      try:
        packfiles.append(file(packfile, 'r'))
      except IOError:
        pass
  else:
    for testpack in build.TestPacks():
      try:
        packfile = '../../test/'+testpack+'TestPack.txt'
        packfiles.append(file(packfile, 'r'))
      except IOError:
        pass
  for packfile in packfiles:
    try:
      for testfile in map(lambda s: s.strip(), packfile.readlines()):
        # Ignore blank lines and repeated tests.
        if testfile and not testfile in testfiles:
          testfiles.add(testfile)
      packfile.close()
    except IOError:
      pass


#print test_cpppath, testsource
#print files, testfiles, testsource

# Add test folders to CPPPATH only for this component
if extra_cpppath:
    env = env.Copy()
    env.Prepend(CPPPATH=extra_cpppath)

# Libraries to link against
all_libs = ['test'+project_name] + chaste_libs + other_libs

# Build the library for this project
project_lib = env.StaticLibrary(project_name, files)

# Build the test library for this project
test_lib = env.Library('test'+project_name, testsource)


# Make test output depend on shared libraries, so if implementation changes
# then tests are re-run.
lib_deps = [project_lib, test_lib] # only this project's libraries
#lib_deps.extend(map(lambda lib: '#lib/lib%s.so' % lib, chaste_libs)) # all Chaste libs

# Collect a list of test log files to use as dependencies for the test
# summary generation
test_log_files = []

# Build and run tests of this component
for testfile in testfiles:
    prefix = testfile[:-4]
    runner_cpp = env.Test(prefix+'Runner.cpp', 'test/' + testfile) 
    env.Program(prefix+'Runner', runner_cpp,
                LIBS = all_libs,
                LIBPATH = ['#/lib', '.'] + other_libpaths)
    if not compile_only:
        log_file = env.File(prefix+'.log')
        #env.Depends(log_file, lib_deps)
        test_log_files.append(log_file)
        env.RunTests(log_file, prefix+'Runner')

Return("test_log_files")
