//
//  Deep_pi_pi_MC.cpp
//  
//
//  Created by Hyde, Charles E. on 5/15/2021.
//
/** @file src/Deep_pi_pi_MC.cxx
 *  @brief
*  \f$l N \to l N \pi \pi\f$
*  lepto-production cross sections (electron or muon) based on
*   ``Angular distributions in hard exclusive production of pion pairs'',
*   B. Lehmann-Dronke, A. Schaefer, M. V. Polyakov, and K. Goeke,
*   PHYSICAL REVIEW D, \f$\mathbf{63}\f$ (2001) 114001.
 *
 *  Input beam kinematics allow fixed target, head-on collisions, and crossing-angle collisions.
 *
 *  Created by Hyde, Charles E. on 5/15/2021.
 *
 *   Compile using make Makefile in bash shell (directory src/).
*/

/** @file src/Deep_pi_pi_MC.hpp
 *  @brief File defines global variables (constants, and event-by-event variables).
 *
 */

#include "Deep_pi_pi_MC.hpp"
#include <iostream>

int Init(char *inFILE);

/** @brief
 *  Monte-Carlo Driver
 */
/*
int main(int argc, char *argv[]){
    char FILEin[64];
    double kCM, Q2Max0, yMax0; // Lepton-Ion CM Lepton energy
    sprintf(FILEin,"deep_sigma.inp");
    if (argc>0) {
        printf("Executing %s \n",argv[0]);
    }
    if (argc>1) {
        sprintf(FILEin,"%s",argv[1]);
        printf(" Parameter input file = ");
        for (int loop=1; loop<argc; loop++) {
            printf("%s   ",argv[loop]);
        }
        printf(" \n");
    }
 */
int main(char FILEin[64] = "deep_pi_pi_JLEIC.txt"){
    //char FILEin[64] = "deep_pi_pi_JLEIC.txt";
    if (Init(FILEin)!=1){
        printf("Invalid initialization \n");
        return -1;
    }

    return 1;
} // main

/** @brief
 *  Read input file to initialize Monte-Carlo event generation.
 */
int Init(char *inFile){
    FILE  * inF;
    inF=fopen(inFile,"readonly");
    printf("Input control parameters\n");
    const int nchar=132;
    double vx,vy,vz, ex,ey,ep, bx,by ;
    double ELepton, EIon;
    char line[nchar], lepton[15], ion[15];
    mPion = dbPDG->GetParticle(211)->Mass();
    //  Read header line
    fgets(line,nchar,inF); printf("%s",line);
    //
    //  Read lepton header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read lepton beam descriptor
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %s", &vx, &vy, &vz, &ex, &ey, &ep, &bx, &by, lepton);
    if (strncmp(lepton,"electron",7)==0) {
        mLepton = dbPDG->GetParticle(11)->Mass();
    } else if (strncmp(lepton,"muon",4)==0){
        mLepton = dbPDG->GetParticle(13)->Mass();
    } else {
        printf("invalid lepton = %s \n", lepton);
        return -1;
    }
    long double ELSq = (vx*vx) + (vy*vy) + (vz*vz) + mLepton*mLepton;
    ELepton = sqrt(ELSq);
    printf("(vx,vy,vz) = (%10.4f,%10.4f,%10.4f), mLepton = %10.3e, ELepton = %13.8lf, ELSq=%13.8Lf \n",
           vx,vy,vz, mLepton, ELepton, ELSq);
    k4Beam0.SetPxPyPzE(vx,vy,vz,ELepton);
    emitt_e[0] = ex*mLepton/ELepton; emitt_e[1] = ey*mLepton/ELepton; emitt_e[2] = ep;
    betaIP_e[0] = bx;                betaIP_e[1] = by;
    if (ex+ey+ep <= 0.0){
        eSmear = false;
    }
    //
    //  Read Ion header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read Ion Beam descripton
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf %s", &vx, &vy, &vz, &ex, &ey, &ep, &bx, &by, ion);
    if (strncmp(ion,"proton",6)==0){
        MIon = dbPDG->GetParticle(2212)->Mass();
    } else if (strncmp(ion,"neutron",7)==0) {
        MIon = dbPDG->GetParticle(2112)->Mass();
    } else if (strncmp(ion,"nucleon",7)==0) {
        MIon = ( dbPDG->GetParticle(2212)->Mass() + dbPDG->GetParticle(2112)->Mass() )/2.0;
    } else {
        printf("invalid ion = %s \n", ion);
        return -1;
    }
    EIon = sqrt(MIon*MIon + vx*vx + vy*vy + vz*vz);
    P4Beam0.SetPxPyPzE(vx,vy,vz, EIon);
    //  Calculate Geometrical (Transverse) Emittance
    emitt_i[0] = ex*MIon/EIon; emitt_i[1]=ey*MIon/EIon; emitt_i[2] = ep;
    betaIP_i[0] = bx;                betaIP_i[1] = by;
    if (ex+ey+ep <= 0.0){
        iSmear = false;
    }
    printf("LeptonMassSq, k4Beam0.M2(), P4Beam0.M2() = %13.6e, %13.6e, %10.6f \n",
           mLepton*mLepton, k4Beam0.M2(), P4Beam0.M2());
    
    //  Read Kinematics header line
    fgets(line,nchar,inF); printf("%s",line);
    //  Read MC Kinematics bounds
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf %d", &vx, &vy, &ex, &ey, &nEvents);
    Q2Min = vx;
    Q2Max = vy;
    yMin  = ex;
    yMax  = ey;
    //  Hadronic Final state bounds
    fgets(line,nchar,inF); printf("%s",line);
    fgets(line,nchar,inF); printf("%s",line);
    sscanf(line,"%lf %lf %lf %lf", &vx, &vy, &ex, &ey);
    MpipiMin  = vx;
    MpipiMax  = vy;
    csPiPiMin = ex;
    csPiPiMax = ey;
    if (MpipiMin<2.*mPion) {
        MpipiMin = 2.0*mPion;
    }
    if (MpipiMax<MpipiMin) {
        return -1;
    }
    if (csPiPiMax>1.0){
        csPiPiMax=1.0;
    }
    if (csPiPiMin<-1.0){
        csPiPiMin=-1.0;
    }
    

    printf("Generate %d MC events \n", nEvents);
    //  Setup up nominal beam light-cone vectors and transverse coordinates
    k4Beam = k4Beam0;
    P4Beam = P4Beam0;
//    EventLightCone();
    n4_e0      = n4_e;
    n4Tilde_e0 = n4Tilde_e;
    X4_e0      = X4_e;
    Y4_e0      = Y4_e;
    W2Threshold= pow(MIon+2.0*mPion,2);
    fclose(inF);
    return 1;
}   //Init()


