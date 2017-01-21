import os
import sys
import pulsar as psr

def run(mm):
    tester=psr.PyTester("Testing the Pulsar-ForceManII API")
    wfn=psr.Wavefunction()
    wfn.system=psr.make_system("""
    0 1
    O     -2.217302    2.805476    2.146856
    H     -2.260473    3.658861    2.578259
    H     -2.256339    3.006951    1.211915
    """)
    atom_types=[2001,2002,2002]
    ff="AMBER99"
    model="HARMONIC_OSCILLATOR"
    coord="BOND"
    mod_key="H.O."
    mm.load_module("fmanii","FFTerm",mod_key)
    mm.change_option(mod_key,"FORCE_FIELD",ff)
    mm.change_option(mod_key,"ATOM_TYPES",atom_types)
    mm.change_option(mod_key,"MODEL_NAME",model)
    mm.change_option(mod_key,"COORD_NAME",coord)
    my_mod=mm.get_module(mod_key,0)
    newwfn,deriv=my_mod.deriv(0,wfn)
    print(deriv)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)

