# -*-Makefile-*-

# Needed to make `"$@"` usable in recipes
set positional-arguments := true

default:
  just run --beam-on 10

setup BUILD SRC:
  #!/usr/bin/env sh
  if [[ ! -d {{BUILD}} ]]; then
    meson setup {{BUILD}} {{SRC}}
  fi

install-tests:
  just install-lib-quiet
  just setup build/neutron_bg-tests test
  meson compile -C build/neutron_bg-tests
  meson install -C build/neutron_bg-tests

test *REST: install-tests
  sh install/neutron_bg-tests/run-each-test-in-separate-process.sh "$@"

catch2-demo *REST:
  echo "$@"
  meson setup      build/neutron_bg-tests test
  meson compile -C build/neutron_bg-tests
  meson install -C build/neutron_bg-tests
  install/neutron_bg-tests/bin/catch2-demo-test "$@"

build:
  just setup build/neutron_bg src
  meson compile -C build/neutron_bg

install: build
  meson install -C build/neutron_bg


install-lib-quiet:
  #!/usr/bin/env sh
  if [[ -f install/neutron_bg/lib/libneutron_bg.so    ||
        -f install/neutron_bg/lib/libneutron_bg.dylib ]]; then
    install_output=$(just install)
    install_ok=$?
    if ! test install_ok; then
      echo $install_output
      exit $install_ok
    fi
  else
    just install
  fi

run *ARGS: install
  #!/usr/bin/env sh
  sh execute-with-nixgl-if-needed.sh ./install/neutron_bg/bin/neutron_bg "$@"
  exit $?

clean:
  rm build install -rf
