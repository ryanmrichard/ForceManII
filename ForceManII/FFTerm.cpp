#include "ForceManII/FFTerm.hpp"
#include "ForceManII/Common.hpp"

using namespace std;
namespace FManII {

const Vector d1(const Vector& Carts,
                const vector<IVector>& ans,
                const InternalCoordinates& coord,
                const Vector& dm)
{
    DEBUG_CHECK(ans.size()==dm.size(),"Derivative sizes are incompatible");
    Vector deriv(Carts.size());
    for(size_t coordi=0;coordi<ans.size();++coordi)
    {
        const IVector atoms=ans[coordi];
        const Vector dc=coord.deriv(1,Carts,atoms);
        for(size_t i=0;i<atoms.size();++i)
            for(size_t j=0;j<3;++j)
                deriv[atoms[i]*3+j]+=dm[coordi]*dc[i*3+j];
    }
    return deriv;
}

Vector FFTerm::deriv(size_t order,const map<string,Vector>& ps,
                          const Molecule& cs)const{
    CHECK(order<2,"Higher order derivatives are not yet implemented!!!");
    vector<Vector> incoords;
    incoords.push_back(cs.coords.at(coord_->name));
    const Vector dm=model_->deriv(order,ps,incoords);
    if(order==0)return dm;
    const Vector dc1=d1(*cs.carts,cs.atom_numbers.at(coord_->name),*coord_,dm);
    if(order==1)return dc1;
}
}
