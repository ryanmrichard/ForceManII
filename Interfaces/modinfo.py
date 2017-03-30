from pulsar import OptionType

minfo={
"ForceField":{
"type":"c_module",
"base":"EnergyMethod",
"modpath":"fmanii_pulsar.so",
"version":"0.1",
"authors":"Ryan M. Richard",
"refs":"",
"description":"Runs an entire FF and returns the deriv",
"options":
{
    "FORCE_FIELD" : (OptionType.String,None,True,None,'The ff to use.'),
    "ATOM_TYPES"  : (OptionType.ListInt,None,True,None,"The atom types"),
    "SKIP_MODELS" : (OptionType.ListString,[],False,None,
        "Optional list of models not to include"),
    "SKIP_COORDS" : (OptionType.ListString,[],False,None,
        "Optional list of coords not to include"),
    "MAX_DERIV"   : (OptionType.Int,1,False,None,
        "The maximum analytic derivative")
}
},

"FFTerm":
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
        "MODEL_NAME"  : (OptionType.String,None,False,None,"The model"),
        "COORD_NAME"  : (OptionType.String,None,False,None,"The coord type"),
        "MAX_DERIV"   : (OptionType.Int,1,False,None,"The maximum analytic derivative")
    }
},
"FFCharges":{
    "type":"c_module",
    "base":"PropertyCalculator",
    "modpath":"fmanii_pulsar.so",
    "version":"0.1",
    "authors":"Ryan M. Richard",
    "refs":"",
    "description":"Returns the charge of each atom",
    "options":
    {
        "FORCE_FIELD" : (OptionType.String,None,True,None,'The ff to use.'),
        "ATOM_TYPES"  : (OptionType.ListInt,None,True,None,"The atom types"),
    }
},
}
