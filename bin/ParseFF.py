""" Given a force field in Tinker format this script will write out a C++ source
file for that force field
"""
from FManII import *

class ForceField:
    def __init__(self):
        self.rad_type=""
        self.rad_size=""
        self.params={}
        self.type2class={}
        self.terms={}
        self.orderrules={}
        self.paramtypes={}
        self.combrules={}
        self.scale_factors={}

def check_add(p,t,pt,k,v):
    if t not in p:p[t]={}
    if pt not in p[t]:p[t][pt]={}
    p[t][pt][k]=v

def read_ff(ff_file):
    """Parses a tinker style force field file"""
    ff=ForceField()
    with open(ff_file,"r") as f:
        for line in f:
            da_line=line.split()
            if len(da_line)<1:continue
            if da_line[0]=="radiusrule":
                ff.combrules[(models["lj"],params["sigma"])]=comb_rules[da_line[1]]
            elif da_line[0]=="radiustype":ff.rad_type=da_line[1]
            elif da_line[0]=="radiussize":ff.rad_size=da_line[1]
            elif da_line[0]=="torsionunit":
                ff.scale_factors[ffterms["ft"]]=da_line[1]
            elif da_line[0]=="imptorunit":
                ff.scale_factors[ffterms["fi"]]=da_line[1]
            elif da_line[0]=="vdwindex":
                idxtype=typetypes[da_line[1].lower()]
                ff.paramtypes[ffterms["lj14"]]=idxtype
                ff.paramtypes[ffterms["lj"]]=idxtype
            elif da_line[0]=="epsilonrule":
                ff.combrules[(models["lj"],params["epsilon"])]=comb_rules[da_line[1]]
            elif da_line[0]=="vdw-14-scale":
                ff.scale_factors[(models["lj"],intcoords["pair14"])]=str(1.0/float(da_line[1]))
            elif da_line[0]=="chg-14-scale":
                ff.scale_factors[(models["cl"],intcoords["pair14"])]=str(1.0/float(da_line[1]))
            elif da_line[0]=="electric":ff.electric=da_line[1]
            elif da_line[0]=="dielectric":ff.dielectric=da_line[1]
            elif da_line[0]=="atom":ff.type2class[int(da_line[1])]=int(da_line[2])
            if da_line[0]=="bond":
                i=int(da_line[1]);j=int(da_line[2])
                pair=(i,j) if i<=j else (j,i)
                ffterm=(models["ho"],intcoords["bond"])
                ff.orderrules[ffterm]="pair_order"
                ff.paramtypes[ffterm]=typetypes["class"]
                ff.terms[ffterm]="HarmonicBond"
                #Tinker bakes the 1/2 into k already
                k,r0=2*float(da_line[3])*k2au,float(da_line[4])*ang2au
                check_add(ff.params,ffterm,params["k"],pair,[k])
                check_add(ff.params,ffterm,params["r0"],pair,[r0])
            if da_line[0]=="angle":
                i=int(da_line[1]);j=int(da_line[2]);k=int(da_line[3])
                triple=(i,j,k) if i<=k else (k,j,i)
                ff.orderrules[ffterms["ha"]]="angle_order"
                ff.paramtypes[ffterms["ha"]]=typetypes["class"]
                ff.terms[ffterms["ha"]]="HarmonicAngle"
                k,r0=2*float(da_line[4])*anglek2au,float(da_line[5])*deg2rad
                check_add(ff.params,ffterms["ha"],params["k"],triple,[k])
                check_add(ff.params,ffterms["ha"],params["r0"],triple,[r0])
            if da_line[0]=="torsion":
                i,j,k,l=int(da_line[1]),int(da_line[2]),\
                        int(da_line[3]),int(da_line[4])
                quad=(i,j,k,l) if (j<k or (j==k and i<l)) else (l,k,j,i)
                ffterm=(models["fs"],intcoords["torsion"])
                ff.orderrules[ffterm]="torsion_order"
                ff.paramtypes[ffterm]=typetypes["class"]
                ff.terms[ffterm]="FourierTorsion"
                v,phi,n=[],[],[]
                for ti in range(3):
                    if len(da_line)<(6+3*ti):break
                    v.append(float(da_line[5+3*ti])*kcalmol2au)
                    phi.append(float(da_line[6+3*ti])*deg2rad)
                    n.append(float(da_line[7+3*ti]))
                check_add(ff.params,ffterm,params["v"],quad,v)
                check_add(ff.params,ffterm,params["phi"],quad,phi)
                check_add(ff.params,ffterm,params["n"],quad,n)
            if da_line[0]=="imptors":
                i,j,k,l=int(da_line[1]),int(da_line[2]),int(da_line[3]),int(da_line[4])
                jkp=[i,j,l]
                jkp.sort()
                quad=(jkp[0],k,jkp[1],jkp[2])
                v,phi,n=float(da_line[5])*kcalmol2au,\
                        float(da_line[6])*deg2rad,\
                        float(da_line[7])
                ffterm=(models["fs"],intcoords["imp"])
                ff.paramtypes[ffterm]=typetypes["class"]
                ff.orderrules[ffterm]="imp_order"
                ff.terms[ffterm]="FourierImproperTorsion"
                check_add(ff.params,ffterm,params["v"],quad,[v])
                check_add(ff.params,ffterm,params["phi"],quad,[phi])
                check_add(ff.params,ffterm,params["n"],quad,[n])
            if da_line[0]=="charge":
                i=int(da_line[1])
                elem=(i,)
                ffterm=(models["cl"],intcoords["pair14"])
                ff.paramtypes[ffterm]=typetypes["type"]
                ff.terms[ffterm]="Electrostatics14"
                ffterm=(models["cl"],intcoords["pair"])
                ff.paramtypes[ffterm]=typetypes["type"]
                ff.terms[ffterm]="ElectrostaticsPair"
                check_add(ff.params,ffterm,params["q"],elem,[float(da_line[2])])
                ff.combrules[(models["cl"],params["q"])]=comb_rules["PRODUCT"]
            if da_line[0]=="vdw":
                i=int(da_line[1])
                elem=(i,)
                val=float(da_line[2])*ang2au
                if ff.rad_size=="RADIUS":val*=2.0
                if ff.rad_type=="SIGMA":val*=(2.0**(1/6))
                ep=float(da_line[3])*kcalmol2au
                check_add(ff.params,ffterms["lj"],params["sigma"],elem,[val])
                check_add(ff.params,ffterms["lj"],params["epsilon"],elem,[ep])
                for ptype in ["pair14","pair"]:
                    ffterm=(models["lj"],intcoords[ptype])
                    if ffterm not in ff.paramtypes:
                        ff.paramtypes[ffterm]=typetypes["class"]
                    ff.terms[ffterm]="LJ14" if ptype=="pair14" else "LJPair"
            if da_line[0]=="vdw14":
                i=int(da_line[1])
                elem=(i,)
                val=float(da_line[2])*ang2au
                if ff.rad_size=="RADIUS":val*=2.0
                if ff.rad_type=="SIGMA":val*=(2.0**(1/6))
                ep=float(da_line[3])*kcalmol2au
                check_add(ff.params,ffterms["lj14"],params["sigma"],elem,[val])
                check_add(ff.params,ffterms["lj14"],params["epsilon"],elem,[ep])


    return ff



def print_parms(f,mol_name,parms):
    """Prints parameter map in initializer list format
    Theactual map is coordtype by parameter type by atom types"""
    var_name=mol_name.lower()+"_FF_params"
    f.write("{\n")#outer most map
    for ct,rest in parms.items():
        f.write("{{"+ct[0]+","+ct[1]+"},{\n")#Start of pair 1st map and start of middle map
        for pt,rest1 in rest.items():
            f.write("   {"+pt+",{\n")#Start of pair middle map and start of last map
            for atoms,value in rest1.items():
                f.write("      {{")#Start of pair last map and start of atom vector
                for ai in atoms:f.write(str(ai)+",")
                f.write("},{")#End atom vector start param vector
                for pi in value:f.write(str(pi)+",")
                f.write("}},\n")#end param vector, end pair last map
            f.write("      }},\n")#end last map end pair middle map
        f.write("   }},\n")#end middle map end 1st pair
    f.write("},\n")#end outer map

##################################### Actual Script ############################
def main():
    if len(sys.argv) !=2:
        raise RuntimeError("Usage: python3 ParseFF.py <param_file>")
    ff_file=sys.argv[1];
    if not os.path.isfile(ff_file):
        raise RuntimeError("Parameter file does not exist")
    ff=read_ff(ff_file)
    ff_name=os.path.splitext(ff_file)[0]

    f=open(ff_name+".cpp","w")
    f.write("//This file is autogenerated from ParseFF.py\n")
    f.write("//If you value your sanity do not try to manually edit it!!!\n")
    f.write("//Owing to GCC's trouble with initializer lists we avoid using\n")
    f.write("//them to initialize our const objects\n\n")
    f.write("#include<ForceManII/FManII.hpp>\n")
    f.write("static FManII::ForceField make_ff(){\n")
    f.write("FManII::ForceField ff;\n")
    new_terms={}
    for i,key in enumerate(ff.terms.items()):
        new_terms[key[0]]="___FFTERM"+str(i)+"___"
        f.write("const FManII::FFTerm_t "+new_terms[key[0]]+"=std::make_pair("+key[0][0]+","+key[0][1]+");\n")
    for i,it in enumerate(params.items()):
        new_terms[it[1]]="___PARAM"+str(i)+"___"
        f.write("const std::string& "+new_terms[it[1]]+"="+it[1]+";\n")
    for ct,rest in ff.params.items():
        for pt,rest1 in rest.items():
            for atoms,value in rest1.items():
                f.write("ff.params.add_param("+new_terms[ct]+","+new_terms[pt]+",")
                if(len(atoms)>1):
                    f.write("std::vector<size_t>({")
                    for ai in atoms:f.write(str(ai)+",")
                    f.write("}),")
                else: f.write(str(atoms[0])+",")
                if(len(value)>1):
                    f.write("std::vector<double>({")
                    for pi in value:f.write(str(pi)+",")
                    f.write("})")
                else:f.write(str(value[0]))
                f.write(");\n")
    for key,value in ff.type2class.items():
        f.write("ff.type2class["+str(key)+"]="+str(value)+";\n")
    for key,value in ff.terms.items():
        f.write("ff.terms.emplace("+new_terms[key]+",std::move(FManII::"+value+"()));\n")
    for key,value in ff.orderrules.items():
        f.write("ff.orderrules["+new_terms[key]+"]=FManII::"+value+";\n")
    for key,value in ff.paramtypes.items():
        f.write("ff.paramtypes["+new_terms[key]+"]="+value+";\n")
    for key,value in ff.combrules.items():
        f.write("ff.combrules[std::make_pair("+key[0]+","+key[1]+")]="+value+";\n")
    for key,value in ff.scale_factors.items():
        f.write("ff.scale_factors["+new_terms[key]+"]="+value+";\n")
    f.write("ff.link_terms("+new_terms[ffterms["lj14"]]+
            ","+new_terms[ffterms["lj"]]+");\n")
    f.write("ff.link_terms("+new_terms[ffterms["cl14"]]+","+
            new_terms[ffterms["cl"]]+");\n")
    f.write("return ff;\n}\n")
    f.write("const FManII::ForceField FManII::"+ff_name+"=make_ff();\n")
    f.close()
if __name__ == "__main__":
    main()
