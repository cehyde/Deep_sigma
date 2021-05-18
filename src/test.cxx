//  Testing the use of TLorentzVector as argument of a function call

#include <TLorentzVector.h>

void Civita(TLorentzVector v1, TLorentzVector v2, TLorentzVector v3, double *out4[4] );

int test(){
    TLorentzVector T4, X4, Y4, Z4
    X4.SetXYZT(1.0,0.0,0.0,0.0);
    Y4.SetXYZT(0.0,1.0,0.0,0.0);
    Z4.SetXYZT(0.0,0.0,1.0,0.0);
    double get4[4];
    Civita4(X4,Y4,Z4,get4);
    
    return 1;
    
}
