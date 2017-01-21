from pulsar import OptionType

minfo={"FFTerm":
{
    "type":"c_module",
    "base":"EnergyMethod",
    "modpath":"fmanii_pulsar.so",
    "version":"0.1",
    "authors":"Ryan M. Richard",
    "refs":"",
    "description":"Runs a single FF term and returns the deriv",
    "options":
    {
        "FORCE_FIELD" : (OptionType.String,None,True,None,'The ff to use.'),
        "ATOM_TYPES"  : (OptionType.ListInt,None,True,None,"The atom types"),
        "MODEL_NAME"  : (OptionType.String,None,True,None,"The model"),
        "COORD_NAME"  : (OptionType.String,None,True,None,"The coord type"),
    }
}}
