project( 'neutron_bg-tests'
       , 'cpp'
       , version : 'v0.0.0'
       , default_options : [ 'buildtype=debugoptimized'
                           , 'debug=true'
                           , 'optimization=2'
                           , 'cpp_std=c++20'
                           , 'prefix=@0@/../install/neutron_bg-tests/'.format(meson.source_root())
                           ]
       )

install_prefix = get_option('prefix')

neutron_bg_lib  = dependency( 'neutron_bg'
                   , method  : 'pkg-config'
                   , required: true
                   )

nain4  = dependency( 'nain4'
                   , method  : 'pkg-config'
                   , required: true
                   )

geant4 = dependency( 'Geant4'
                   , method  : 'cmake'
                   , required: true
                   , components: ['ui_all', 'vis_all']
                   , modules    : nain4.get_variable(pkgconfig: 'geant4_modules').split(','))

catch2 = dependency( 'catch2'
                   , method   : 'pkg-config'
                   , required : true)

test_deps       = [neutron_bg_lib, nain4, geant4, catch2]
test_sources    = ['catch2-main-test.cc', 'test-LXe.cc']
c2_demo_sources = ['catch2-main-test.cc', 'test-catch2-demo.cc']

geant4_include = geant4.get_variable(cmake    : 'Geant4_INCLUDE_DIRS')
nain4_include  =  nain4.get_variable(pkgconfig: 'includedir'         )
neutron_bg_lib_include  = neutron_bg_lib.get_variable(pkgconfig: 'includedir'         )

test_executable = 'neutron_bg-test'
demo_executable = 'catch2-demo-test'

test_exe = executable( test_executable
                     , test_sources
                     , include_directories: [neutron_bg_lib_include, nain4_include, geant4_include]
                     , dependencies       : test_deps
                     , install            : true
                     )

c2_demo_exe = executable( demo_executable
                        , c2_demo_sources
                        , include_directories: []
                        , dependencies       : [catch2]
                        , install            : true
                        )

conf_data = configuration_data({
'install_prefix'  : install_prefix,
'test_executable' : test_executable,
})

configure_file(
  input         : '../run-each-test-in-separate-process.sh.in',
  output        :    'run-each-test-in-separate-process.sh',
  configuration : conf_data,
  install       : true,
  install_dir   : install_prefix
)
