{pkgs} :
{
  shell-shared = {
    G4_DIR = "${pkgs.geant4}";

    shellHook = ''
      export neutron_bg_lib_PATH=$PWD/install/neutron_bg/lib
      export LD_LIBRARY_PATH=$neutron_bg_lib_PATH:$LD_LIBRARY_PATH;
      export PKG_CONFIG_PATH=$neutron_bg_lib_PATH/pkgconfig:$PKG_CONFIG_PATH;
    '';
  };
}
