//
//  Deep_event.cxx
//  
//
//  Created by Charles Earl Hyde on 12/4/17.
//  Revised 16 May 2021 to replace qft++ structures
//
//
/** @file include/Deep_event.cxx
 *  @brief File contains event-by-event generation functions for  \f$e p \to e' p' \pi \pi\f$
 */



#ifndef Deep_event_h
#define Deep_event_h

/** @brief
 * Initialize lightcone vectors n4_e, n4Tilde_e, X4_e, Y4_e.
 *
 *  \f{align}{
 *  n_e^\mu &=  \displaystyle \left[k^\mu\left(1+\sqrt{1-\delta_l}\right) - \frac{m_l^2}{k\cdot P} P^\mu \right]
 *    \left/\left[2 \sqrt{(k\cdot P)(1-\delta_l)} \right]\right.,
 *    \qquad\qquad k^2 = m_l^2, \qquad \delta_l = \left[\frac{m_l M_{\rm Ion}}{k\cdot P}\right]^2 \\
 *  \widetilde n_e^\mu &= \left[P^\mu 
 *                    - \frac{M_{\rm Ion}^2}{(k\cdot P)\left(1+\sqrt{1-\delta_l}\right)}k^\mu
 *                     \right]\left/\frac{1}{\sqrt{(k\cdot P)(1-\delta_l)}} \right. \\
 *  n_e\cdot n_e &= 0 = \widetilde n_e\cdot \widetilde n_e \\
 *  n_e\cdot \widetilde n_e &= 1 \\
 *  Transverse vectors \\
 *  {\tt Y4Det}^\mu  = [0,0,1,0] = up in Detector frame
 *  [X4_{e0}]_\sigma   =\epsilon_{\mu\nu\rho\sigma}n4_e^\mu n4Tilde_e^\nu Y4_Det^\rho \\
 *  X4_e^\mu         = X4_{e0}^\mu/\sqrt{-X4_{e0}\cdot X4_{e0} } \\
 *  [Y4_e]_\sigma      = \epsilon_{\mu\nu\rho\sigma}n4_e^\mu X4_e^\nu n4Tilde_e^\rho\\
 * \f}
 * If beam emmitance values in input file are positive, incident beam 4-vectors
 * \f$k^\mu, P^\mu\f$
 * are generated with gaussian longitudinal and transverse emmittance relative to nominal input values
 * \f$k_0^\mu, P_0^\mu\f$.
 *
 * The transverse 4-vectors \f$ X_e^\mu,\, Y_e^\mu\f$ are defined assuming
 * neither incident beam can ever be in the vertical direction.
 */
int LeviCivita4vec(double vec1[4], double vec2[4], double vec3[4], double *vec4out)
{
    /** @brief
           * Construct a 4-vector contraction
     *\f{align}{
     *{\tt vec4out}_\mu &= \epsilon_{\mu \nu \rho \sigma} {\tt vec1}^\nu {\tt vec2}^\rho {\tt vec3}^\sigma\\
     * \epsilon_{0123} &= 1\\
     * \f}
     */
    int ndim = 4;
    int yes  = 0;
    double civita = 1.0;
    for (int ii=0; ii<ndim; ii++) {
        vec4out[ii] = 0.0;
        for (int jj=0; jj<ndim; jj++) {
            if (jj!=ii){
                for (int kk=0; kk<ndim; kk++) {
                    if((kk!=jj)&&(kk!=ii)){
                        for (int ll=0; ll<ndim; ll++) {
                            if((ll!=kk)&&(ll!=jj)&&(ll!=ii)){
                                vec4out[ii] += civita*vec1[jj]*vec2[kk]*vec3[ll];
                                civita*=-1.0;
                            } // ll!= {kk, jj, ii}
                        } // ll for loop
                        civita*=-1.0;
                    } //kk!={jj, ii}
                } // kk for loop
                civita *= -1.0;
            } //jj!=ii
        } // jj for loop
        civita *= -1.0;
    } // ii for loop
    yes = 1;
    return yes;
} // LeviCivita4Vec

//double LeviCivitaScalar(double *vec1, double *vec2, double *vec3, double &vec4);

double LeviCivitaScalar(TLorentzVector vec4_1, TLorentzVector vec4_2, TLorentzVector vec4_3, TLorentzVector vec4_4){
/** @brief
       * Construct a the scaler anti-symmetric contraction of four space-time vectors
 *\f{align}{
 *{\tt Scalar} &= \epsilon_{\mu \nu \rho \sigma} {\tt vec4_1}^\mu {\tt vec4_2}^\nu {\tt vec4_3}^\rho {\tt vec4_4}^\sigma\\
 * \epsilon_{0123} &= 1\\
 * \f}
 */
    double Scalar=0.0;
    int ndim = 4;
    double term;
    double civita = 1.0;
    for (int ii=0; ii<ndim; ii++) {
        for (int jj=0; jj<ndim; jj++) {
            if (jj!=ii){
                for (int kk=0; kk<ndim; kk++) {
                    if((kk!=jj)&&(kk!=ii)){
                        for (int ll=0; ll<ndim; ll++) {
                            if((ll!=kk)&&(ll!=jj)&&(ll!=ii)){
                                Scalar += civita*vec4_1[ii]*vec4_2[jj]*vec4_3[kk]*vec4_4[ll];
                                civita*=-1.0;
                            } // ll!= {kk, jj, ii}
                        } // ll for loop
                        civita*=-1.0;
                    } //kk!={jj, ii}
                } // kk for loop
                civita *= -1.0;
            } //jj!=ii
        } // jj for loop
        civita *= -1.0;
    } // ii for loop
    return Scalar;
}

int EventLightCone(){
    //  Initialize lightcone vectors n4, n4Tilde
    //  n4 = (k4Beam - beta*PBeam)/norm1
    //  n4Tilde = (PBeam+betaPr*k4Beam)/norm2;
    //  n4\cdot n4 = 0 = n4Tilde\cdot n4Tilde
    //  n4\cdot n4Tilde = 1;
    sqrtDL      = mLepton*MIon/(k4Beam.Dot(P4Beam));
    deltaL      = sqrtDL*sqrtDL;

    double sqrt_One_d  = sqrt(1.0-deltaL);
    double norm        = 1.0/(sqrt(k4Beam.Dot(P4Beam))*sqrt_One_d);
    TLorentzVector   temp4, X4,Y4,Z4,T4;
    X4.SetPxPyPzE(1.0,0.0,0.0,0.0);
    Y4.SetPxPyPzE(0.0,1.0,0.0,0.0);
    Z4.SetPxPyPzE(0.0,0.0,1.0,0.0);
    T4.SetPxPyPzE(0.0,0.0,0.0,1.0);
    printf("Test of LeviCivita Contraction = %13.6f  \n",
           LeviCivitaScalar(T4,X4,Y4,Z4));
    static  bool first = true;
    static int nSpaceTime=4;
    double vec1[nSpaceTime], vec2[nSpaceTime], vec3[nSpaceTime], vec4out[nSpaceTime];
    n4_e       = k4Beam;
    n4_e      *= (1.0+sqrt_One_d);
    temp4      = P4Beam;
    temp4     *= -sqrtDL;
    n4_e      += temp4;
    n4_e      *= norm/2.0; //** n4_e \f$= n_e^\mu \f$
    n4Tilde_e  = P4Beam;
    temp4      = k4Beam;
    temp4     *=(-MIon*MIon/((1.+sqrt_One_d)*(k4Beam.Dot(P4Beam))));
    n4Tilde_e += temp4;
    n4Tilde_e *= norm;
    vec1[0] = n4_e.E();
    vec1[1] = n4_e.Px();
    vec1[2] = n4_e.Py();
    vec1[3] = n4_e.Pz();
    vec2[0] = n4Tilde_e.E();
    vec2[1] = n4Tilde_e.Px();
    vec2[2] = n4Tilde_e.Py();
    vec2[3] = n4Tilde_e.Pz();
    vec3[0] = Y4_Det.E();
    vec3[1] = Y4_Det.Px();
    vec3[2] = Y4_Det.Py();
    vec3[3] = Y4_Det.Pz();
    // Check sign of Y4_e and X4_e!!
    LeviCivita4vec(vec1,vec2,vec3,vec4out);
    X4_e.SetPxPyPzE(-vec4out[1],-vec4out[2],-vec4out[3],vec4out[0]);
    norm = X4_e.M2();
    if (norm>=0.0) {
        printf("Error space-like X4_e.M2() = %8.3f \n", norm);
        return -1;
    }
    X4_e      *= 1./sqrt(-norm);
    vec3[0] = X4_e.E();
    vec3[1] = X4_e.Px();
    vec3[2] = X4_e.Py();
    vec3[3] = X4_e.Pz();
    LeviCivita4vec(vec1,vec3,vec2,vec4out);
    Y4_e.SetPxPyPzE(-vec4out[1],-vec4out[2],-vec4out[3],vec4out[0]);
    
    
    if (first) {
        printf("Lepton mass squared = %10.3g GeV2 \n", k4Beam.M2());
        printf("n4_e^2 = %10.3g = ?0? = %10.3g = n4Tilde_e^2 \n", n4_e.M2(), n4Tilde_e.M2() );
        printf("n4_e * n4Tilde_e = %10.7f =? 1 \n", n4_e.Dot(n4Tilde_e));
        /*
        n4_e.Print();      printf(" = n4_e^{0,1,2,3}      \n");
        n4Tilde_e.Print(); printf(" = n4Tilde_e^{0,1,2,3} \n");
        X4_e.Print();      printf(" = X4_e^{0,1,2,3},     X4_e*X4_e = %8.5f \n", X4_e*X4_e);
        Y4_e.Print();      printf(" = Y4_e^{0,1,2,3},     Y4_e*Y4_e = %8.5f \n", Y4_e*Y4_e);
         */
        first = false;
    }
    return 1;
}

/**  @brief Get_Event() generates ep->e p pi pi events uniformly in phase space
*  \f[\left\{Q^2, y=q\cdot P/(k\cdot P), \phi_e, M_{\pi\pi}^2, \cos(\theta_{\pi\pi}^{CM}), \phi_{\pi\pi}^{CM},
      * \cos\theta_{\pi^+}^{Rest}, \phi_{\pi^+}^{Rest} \right\}
      * \f]
 *  Only basis 4-vectors e.g. \f$ n_q, \widetilde n_q, X_q, Y_q\f$ are boosted,
 *  all other variables are invariants.
*/
int Get_Event(int iEvt){

    // Generate incident electron
    double  pperp, k3Pr_perp, denom;
    double  rn1, fi, norm;
    TLorentzVector  tempX, tempY, tempN, tempNT;
    TLorentzVector P4pi1, P4pi2;
    double EpipiCM, PpipiCM, EPprCM, PPprCM,  determinant, PPr_dot_n, PPr_dot_nT;
    double kPr_dot_n, kPr_dot_nT, kPr_perpSq;
//    TLorentzVector n4_q, n4Tilde_q, X4_q, Y4_q, X4_pipi, Y4_pipi;
    double px, py, pz, pSmear, n0, nx,ny,nz;
    double EIon, ELepton;
    double csPiPiCM, phiPiPiCM;
    TVector3 boostCM;
    static bool first = true;
    
//      double ETot, betaCMx, betaCMy, betaCMz;
    k4Beam = k4Beam0;
    P4Beam = P4Beam0;
    if(first){
        EventLightCone();
        n4_e0 = n4_e;
        n4Tilde_e0 = n4Tilde_e;
        X4_e0 = X4_e;
        Y4_e0 = Y4_e;
        first = false;
    }
    if(iEvt==0) printf("eSmear, iSmear = %d, %d, k4Beam0.M2()=%f \n", eSmear, iSmear, k4Beam0.M2());
    if (eSmear) {
        pSmear = ran3.Gaus(1.0,emitt_e[2]);
        k4Beam *= pSmear;
        rn1    = ran3.Rndm();
        fi     = ran3.Rndm()*TwoPi;
        pperp  = sqrt(-2.0*log(rn1))*k4Beam.P();
        tempX  = X4_e0;
        tempX *= (sqrt(emitt_e[0]/betaIP_e[0])*cos(fi)*pperp);
        tempY  = Y4_e0;
        tempY *= (sqrt(emitt_i[1]/betaIP_e[1])*sin(fi)*pperp);
        k4Beam += tempX + tempY;
        px     = k4Beam.Px();
        py     = k4Beam.Py();
        pz     = k4Beam.Pz();
        ELepton= sqrt(px*px+py*py+pz*pz+mLepton*mLepton);
        k4Beam.SetPxPyPzE(px,py,pz,ELepton);
        // printf("Longitudinal smear factor(e) = %10.6f \n", pSmear);
    }
    P4Beam = P4Beam0;
    if (iSmear) {
        pSmear = ran3.Gaus(1.0,emitt_i[2]);
        P4Beam *= pSmear;
        rn1    = ran3.Rndm();
        fi     = ran3.Rndm()*TwoPi;
        pperp  = sqrt(-2.0*log(rn1))*P4Beam.P();
        /*  Diagnostic, solved by initializing betaIP_i
        printf("rn1 = %f, pperp_x, pperp_y = %f, %f \n", rn1,
               sqrt(emitt_i[0]/betaIP_i[0])*cos(fi)*pperp, sqrt(emitt_i[1]/betaIP_i[1])*sin(fi)*pperp);
         */
        tempX  = X4_e0;
        tempX *= sqrt(emitt_i[0]/betaIP_i[0])*cos(fi)*pperp;
        tempY  = Y4_e0;
        tempY *= sqrt(emitt_i[1]/betaIP_i[1])*sin(fi)*pperp;
        P4Beam += tempX + tempY;
        px     = P4Beam.Px();
        py     = P4Beam.Py();
        pz     = P4Beam.Pz();
        EIon   = sqrt(px*px+py*py+pz*pz+MIon*MIon);
        P4Beam.SetPxPyPzE(px,py,pz,EIon);
    }
    if(iEvt==0) {
        printf("eSmear, iSmear = %d, %d, k4Beam0.M2()=%10.3g \n", eSmear, iSmear, k4Beam0.M2());
        printf("k4Beam.M2()=%10.3g, P4Beam.M2()=%11.7f \n", k4Beam.M2(), P4Beam.M2());
    }

    /**
     * After smearing the incident beam momenta,
     * re-define the lepton-Ion light-cone vectors n4_e, n4Tilde_e, X4_e, Y4_e
     */
    if (EventLightCone()<0) {
        printf("Event %d:  Invalid light cone vectors \n",iEvt);
        return -1;
    }
    //  Generate Electron Scattering Kinematics
    psf      = (Q2Max-Q2Min)*(yMax - yMin)*TwoPi;
    k_dot_P  = k4Beam.Dot(P4Beam);
    sMinusM2 = 2.0*k_dot_P;
    s_e      = sMinusM2 + MIon*MIon + mLepton*mLepton;
    Q2       = Q2Min + (Q2Max-Q2Min)*ran3.Rndm();
    yInv     = yMin  + (yMax - yMin)*ran3.Rndm();
    phi_e    = TwoPi*ran3.Rndm();
    xBj      = Q2/(yInv*sMinusM2);
    W2       = MIon*MIon+ yInv*sMinusM2-Q2;
    //  Check for valid kinematics above \pi\pi threshold
    if (W2<W2Threshold) {
        return -1;
    }
    /**
     *Scattered lepton four-vector k4Scat = \f$k'\f$.
     *
     * \f{align}{
     * \left[\begin{array}{r} Q^2 + 2 m_l^2 \\ 2(k\cdot P)(1-y) \end{array}\right]
     *  &= 2 
     * \left[\begin{array}{rr} k\cdot \widetilde n_e,&\ k\cdot n_e \\ 
     *                         P\cdot \widetilde n_e,& P\cdot n_e \end{array}\right]
     * \left[\begin{array}{r} k' \cdot n_e \\ k' \cdot \widetilde n_e \end{array}\right]
     * \\
     * \left[\begin{array}{r} k' \cdot n_e \\ k' \cdot \widetilde n_e \end{array}\right]
     * &= \frac{1}{(k\cdot \widetilde n_e)(P\cdot n_e) - (k\cdot n_e)(P\cdot \widetilde n_e)}
     * \left[\begin{array}{rr} P\cdot n_e,&\ -k\cdot n_e \\
     *                         -P\cdot \widetilde n_e,& k\cdot \widetilde n_e \end{array}\right]
     * \left[\begin{array}{r} Q^2/2 +  m_l^2 \\ (k\cdot P)(1-y) \end{array}\right]
     * \\
     * \mathbf k_\perp^{\prime 2} &= 2(k'\cdot n_e)(k'\cdot \widetilde n_e) - m_l^2
     * \f}
     */
    deltaQ     = Q2*(MIon/(k_dot_P*yInv))*(MIon/(k_dot_P*yInv));
    determinant= (k4Beam.Dot(n4Tilde_e))*(P4Beam.Dot(n4_e))
                 -(k4Beam.Dot(n4_e))*(P4Beam.Dot(n4Tilde_e));
    kPr_dot_n  = (Q2 + 2.0*mLepton*mLepton)*(P4Beam.Dot(n4_e))/2.0-
                 (k4Beam.Dot(n4_e))*(k4Beam.Dot(P4Beam))*(1.0-yInv);
    kPr_dot_n *= 1./determinant;
    kPr_dot_nT = -(P4Beam.Dot(n4Tilde_e))*(Q2 + 2.0*mLepton*mLepton)/2.0 +
                  (k4Beam.Dot(n4Tilde_e))* k_dot_P *(1.0-yInv);
    kPr_dot_nT*= 1./determinant;
    if ((iEvt%1000)==0) {
        printf("Event %d: (e,e') Kinematics (Q2, 2kP, y, x, W2) = (%10.6f, %10.3f, %10.3f, %10.3f, %10.3f) \n",
               iEvt, Q2, sMinusM2, yInv, xBj, W2);
        printf("      k' Determinant = %13.6e \n", determinant);
    }
    kPr_perpSq = mLepton*mLepton -2.0*kPr_dot_n*kPr_dot_nT;   //  = - (3-vector)^2
    if (kPr_perpSq>0.0) {
        printf("Event %d: Invalid kinematics k'_perp^2 (4-vector) = %10.3g \n",
               iEvt, kPr_perpSq);
        return -1;
    } else {
        k3Pr_perp = sqrt(-kPr_perpSq);
    }
    tempN      = n4_e;
    tempN     *= kPr_dot_nT;
    tempNT     = n4Tilde_e;
    tempNT    *= kPr_dot_n;
    tempX      = X4_e;
    tempX     *= k3Pr_perp*cos(phi_e);
    tempY      = Y4_e;
    tempY     *= k3Pr_perp*sin(phi_e);
    k4Scat     = tempN + tempNT + tempX + tempY;
    q4Virt     = k4Beam - k4Scat;
    P4Tot      = P4Beam + q4Virt;

    /**  Define the lightcone vectors of q + P system
     * \f{align}{
     * n_q^\mu &= \frac{1}{\sqrt{2(1+\delta_Q)} } \left[ 
     *            q^\mu \frac{M_{Ion}}{q\cdot P}
     *          + P^\mu \frac{\delta_Q}{M_{Ion}\left(1+\sqrt{1+\delta_Q}\right)}
     * \right], \qquad
     *   \delta_Q = \frac{Q^2 M_\text{Ion}^2}{(q\cdot P)^2} \\
     * \widetilde n_q^\mu &= \frac{1}{\sqrt{2(1+\delta_Q)} }
     *         \left[ P^\mu \frac{\left(1+\sqrt{1+\delta_Q}\right)}{M_\text{Ion}}
     *              - q^\mu \frac{M_\text{Ion}}{q\cdot P} \right]
     * \f}
     */
    denom      = MIon/(sqrt(2.0*(1.0+deltaQ))*(q4Virt.Dot(P4Beam)));
    tempN      = denom*q4Virt;
    denom      = 1.0/(MIon*sqrt(2.0*(1.0+deltaQ)));
    n4_q       = (deltaQ*denom/(1.0+sqrt(1.0+deltaQ)))*P4Beam;
    n4_q      += tempN;
    n4Tilde_q  = (denom*(1.0+sqrt(1.0+deltaQ)))*P4Beam;
    tempN     *= -1.0;
    n4Tilde_q += tempN;
    Y4_q       = -sin(phi_e) * X4_e + cos(phi_e) * Y4_e;
//    X4_q       = (LeviCivita4|n4_q%n4Tilde_q%Y4_q);
    if (iEvt%1000==0){
        printf("      n4_q*n4_q = %13.6e, n4_q*n4Tilde_q =%13.6e, X4_q*X4_q = %10.5f \n",
               n4_q*n4_q,n4_q*n4Tilde_q, X4_q*X4_q);
    }
    
    //  Hadronic Final State
    psf       *= (MpipiMax*MpipiMax-MpipiMin*MpipiMin);
    MpipiSq    = MpipiMin*MpipiMin +
                 (MpipiMax*MpipiMax-MpipiMin*MpipiMin)*ran3.Rndm();
    if (((MIon+sqrt(MpipiSq))*(MIon+sqrt(MpipiSq)))>W2) {
        return -1;
    }
    psf      *= TwoPi*(csPiPiMax-csPiPiMin);
    csPiPiCM  = csPiPiMin + (csPiPiMax-csPiPiMin)*ran3.Rndm();
    // csPiPiCM  = 1.0;
    phiPiPiCM = TwoPi*ran3.Rndm();
    //  Boost basis vectors, rather than physics vectors
    //  Boost basis 4-vectors from Detector Frame to gamma* + p CM frame
    EpipiCM   = (W2-MIon*MIon+MpipiSq)/(2.*sqrt(W2));
    PpipiCM   = sqrt(EpipiCM*EpipiCM-MpipiSq);
    EPprCM    = (W2+MIon*MIon-MpipiSq)/(2.*sqrt(W2));
    PPprCM    = sqrt(EPprCM*EPprCM-MIon*MIon);
    pperp     = PPprCM*sin(acos(csPiPiCM));
    n4_qCM      = n4_q;
    n4Tilde_qCM = n4Tilde_q;
    X4_qCM      = X4_q;
    Y4_qCM      = Y4_q;
    boostCM = (-1.0)*P4Tot.BoostVector();
    n4_qCM.Boost(boostCM);
    n4Tilde_qCM.Boost(boostCM);
    X4_qCM.Boost(boostCM);
    Y4_qCM.Boost(boostCM);
// X4 and Y4 must be pure 3-vectors in q+P CM frame
    if ((X4_qCM.T()*X4_qCM.T()>1.e-12)||(Y4_qCM.T()*Y4_qCM.T()>1.e-12)) {
        //printf("Event %d, (X4_qCM.T(), Y4_qCM.T() = (%10.3g,%10.3g) \n", iEvt,X4_qCM.T(),Y4_qCM.T());
        printf(" Event %d, X4_qCM ", iEvt);
        X4_qCM.Print();
        printf("\n Event %d, Y4_qCM ", iEvt);
        Y4_qCM.Print();
        printf(" \n");
        return -1;
    }
    /**
     * Recoil proton kinematics:
     * - After boosting to the \f$q+P\f$ CM frame, \f$X_q^\mu,\, Y_q^\mu\f$ are pure space-like
     * (vanishing time components) and the space parts of \f$n_q^\mu,\, \widetilde n_q^\mu\f$
     * are anti-colinear.
     * - In this frame, the incident proton is collinear with the unit three vector
     * - \f$ \mathbf v =\left.\left[
     *       \widetilde n_q.Px(), \widetilde n_q.Py(), \widetilde n_q.Pz()
     *   \right]\right/ \widetilde n_q.E()
     *  \f$
     *  - The recoil nucleon and the \f$\pi\pi\f$ system are back-to-back in \f$q+P\f$ CM frame.
     *
     * - Generate uniform distributions in
     *    - cosine of polar-angle: \f$\cos\left(\theta_{\pi\pi}^{CM}\right) \f$; and
     *    -   azimuthal-angle \f$\phi_{\pi\pi}^{CM}\f$
     * -\f[ \left[ \begin{array}{r} E' \\ \left|\mathbf P^\prime\right| \\
     \end{array} \right]^{CM}
     * = \left[ \begin{array}{r} (W^2+M_\text{Ion}^2-M_{\pi\pi}^2) / (2\sqrt{W^2}) \\
     *     \sqrt{\left[E^{\prime\text{CM}}\right]^2 - M_\text{Ion}^2}  \end{array} \right]
     * \f]
     * Construct CM unit 3-vector in incident ion direction
     *\f[ \mathbf v = \left[ \frac{\widetilde n_q.Px()}{\widetilde n_q.E()},\,
     *    \frac{\widetilde n_q.Py()}{\widetilde n_q.E()},\,
     *    \frac{\widetilde n_q.Pz()}{\widetilde n_q.E()}\right]^{CM}     \f]
     * Define a longitudinal four-vector of the recoil proton in the CM frame:
     *\f[ \left[P_L^{\prime \mu}\right]^{CM} = \left[ E^{\prime},
     -\left|\mathbf P'\right| \cos\left(\theta_{\pi\pi}^\text{CM}\right) \mathbf v \right]^\text{CM}    \f]
     * In any frame, the full recoil proton momentum four-vector is defined by:
     * \f[ P^{\prime\, \mu} = \left(P_L^{\prime\,CM}\cdot n_q^{CM} \right) \widetilde n_q^\mu
     *     + \left(P_L^{\prime\,CM}\cdot \widetilde n_q^{CM} \right)  n_q^\mu
     *     - \left|\mathbf P'\right|^\text{CM} \sin\theta_{\pi\pi}^{CM}
     *         \left[ \cos\phi_{\pi\pi}^{CM} X_q^\mu + \sin\phi_{\pi\pi}^{CM} Y_q^\mu
     *        \right]
     * \f]
     * The momentum four-vector of the \f$\pi\pi\f$ system is simply
     *\f[P_{\pi\pi}^{\mu} = q^\mu + P^\mu - P^{\prime\,\mu} \f]
     */

    n0 = n4Tilde_qCM.E();
    nx = n4Tilde_qCM.Px()/n0;
    ny = n4Tilde_qCM.Py()/n0;
    nz = n4Tilde_qCM.Pz()/n0;

    P4Scat.SetPxPyPzE(-PPprCM*csPiPiCM*nx,-PPprCM*csPiPiCM*ny,-PPprCM*csPiPiCM*nz,EPprCM);
    double Mtemp = P4Scat.M2();
    PPr_dot_n   = P4Scat.Dot(n4_qCM);
    PPr_dot_nT  = P4Scat.Dot(n4Tilde_qCM);
    tempN       = PPr_dot_nT*n4_q;
    tempNT      = PPr_dot_n*n4Tilde_q;
    tempX       = (-pperp*cos(phiPiPiCM))*X4_q;
    tempY       = (-pperp*sin(phiPiPiCM))*Y4_q;
    P4Scat        = tempN + tempNT;
    P4Scat       += tempX + tempY;
    P4pipi      = P4Tot - P4Scat;
    Delta4vec   = P4Scat - P4Beam;
    
    if (iEvt%1000 == 0) {
        printf("Event %d: n4_qCM*X4_qCM =%13.6e = %13.6e = n4Tilde_q*X4_q \n",
               iEvt, n4_qCM*X4_qCM, n4Tilde_q*X4_q);
    }
    
    /**
     * Similar method to construct  four vectors of the final state pions
     * -  Generate the \f$ M_{\pi\pi} \rightarrow \pi \pi\f$ kinematics variables
     *    in the rest frame of the \f$ M_{\pi\pi}\f$ system.
     *     - \f$ \cos\theta_\pi^\text{Rest}, \phi_\pi^\text{Rest} \f$ are generated uniformly
     *     - for charge pion, the \f$\pi^+ \f$ is the pion defined by \f$ \theta_\pi^\text{Rest} \f$
     *     - Is the phase space for a neutral pion \f$ 2\pi \f$ or \f$ 4\pi\f$ ?
     */
      double csTh_pi_Rest = 1.0 - ran3.Rndm();
      double phi_pi_Rest  = TwoPi*(0.500 - ran3.Rndm());
      double E_pi_Rest    = sqrt(MpipiSq)/2.0;
    if (E_pi_Rest<mPion){
        printf("Event %d illegal kinematics E_pi_Rest = %13.6e < mPion = %13.6e \n",
               iEvt, E_pi_Rest, mPion);
        return -1;
    }
      double p_pi_Rest    = sqrt(MpipiSq/4.0 - mPion*mPion);
      double beta_Rest    = p_pi_Rest/E_pi_Rest;
      double beta_pi_pi   = P4pipi.P()/P4pipi.E();
      double E_factor     = (beta_Rest*csTh_pi_Rest*beta_pi_pi+1.0)/2.0;
      double p_factor     = (beta_Rest*csTh_pi_Rest/beta_pi_pi+1.0)/2.0;
    /* First compute longitudinal four-vectors of pion 1 & 2, boosted to detector frame
     */
    P4pi1.SetPxPyPzE(P4pipi.Px()*p_factor, P4pipi.Py()*p_factor, P4pipi.Pz()*p_factor,P4pipi.E()*E_factor);
    E_factor = (1.0 - 1.0*beta_Rest*csTh_pi_Rest*beta_pi_pi)/2.0;
    p_factor = (-beta_Rest*csTh_pi_Rest/beta_pi_pi+1.0)/2.0;
    P4pi2.SetPxPyPzE(P4pipi.Px()*p_factor, P4pipi.Py()*p_factor, P4pipi.Pz()*p_factor,P4pipi.E()*E_factor);
    //P4pi2 = P4pipi - Prpi1;

    Y4_pipi = cos(phiPiPiCM)*Y4_q - sin(phiPiPiCM)*X4_q;
    
    
    //X4_pipi = (LeviCivita4|P4pi2%P4pi1%Y4_pipi) /
     //         sqrt( ((P4pi1*P4pi2)*(P4pi1*P4pi2) - (P4pi2*P4pi2)*(P4pi1*P4pi1)) );
    if (iEvt%1000==0){
        printf("Event %d:  (Longitudinal) Y4_pipi*P4pi1 =%13.6e, Y4_pipi*Prpi2 =%13.6e, P4pipi*Y4_pipi =%13.6e \n",
               iEvt, Y4_pipi*P4pi1, Y4_pipi*P4pi2, Y4_pipi*P4pipi);
    }
    pperp   = p_pi_Rest * sin(acos(csTh_pi_Rest));
    if (abs(Y4_pipi*P4pi1)>1.e-6) {
      //  printf("      Y4_pipi*P4pi1 = %13.6e \n", Y4_pipi*P4pi1);
    }
    P4pi1   = P4pi1 +  (pperp*cos(phi_pi_Rest))*X4_pipi +  (pperp*sin(phi_pi_Rest))*Y4_pipi;
    P4pi2   = P4pi2 + (-pperp*cos(phi_pi_Rest))*X4_pipi + (-pperp*sin(phi_pi_Rest))*Y4_pipi;
    
    if (iEvt%1000 == 0) {
        printf("Event %d:  Y4_pipi^2 = %10.6e =? %10.6e = X4_pipi^2, X4_pipi*Y4_pipi = %13.6e \n",
               iEvt, Y4_pipi*Y4_pipi, X4_pipi*X4_pipi, X4_pipi*Y4_pipi );
        printf("            M_pi_pi^2 = %13.6e, P4pi1*P4pi1 + P4pi2*P4pi2 +2*P4pi1*P4pi2 = %13.6e \n",
             MpipiSq, P4pi1*P4pi1 + P4pi2*P4pi2 +2.0*(P4pi1*P4pi2));
    }
    
    // Kinematic checks
    /*
    n0 = n4_qCM.E();
    if (iEvt%1000 == 0) {

        printf("Event %d: 3vector n^2 = %10.5f, 4-vec X4q^2 = %10.f = %10.f = Y4q^2 \n",
               iEvt, nx*nx+ny*ny+nz*nz, X4_qCM.Mass2(), Y4_qCM.Mass2());
        printf("      CM: (nT,n) = [(%10.5f,%10.5f), (%10.5f,%10.5f), (%10.5f,%10.5f), (%10.5f,%10.5f)] \n",
               n4Tilde_qCM.E(),n0,nx,n4_qCM.Px()/n0, ny,n4_qCM.Py()/n0, nz,n4_qCM.Pz()/n0);
        printf("Event %d: m_pipi mass-squared difference = %13.6e GeV2, Recoil Ion MSq = %13.5e GeV2 \n",
               iEvt, P4pipi.Mass2()-MpipiSq, P4Pr.Mass2());
        printf(" 2[n4_q.E()*n4Tilde_q.E()]_CM = %13.6f \n", 2.*(n4_qCM.E())*(n4Tilde_qCM.E()));
        printf("Event %d: (qn)(PnT)+(qnT)(Pn) = %10.5f = qP = %10.5f, q4*q4/Q2 = %10.6f \n",
               iEvt, (n4_q*q4)*(n4Tilde_q*PBeam)+(n4Tilde_q*q4)*(n4_q*PBeam), q4*PBeam, q4.Mass2()/Q2);
        printf("  (q+P).E()_CM = %10.5f = W =  %10.5f = sqrt((q+P)^2) = %10.5f \n",
               ((n4_q*q4)+(n4_q*PBeam))*n4Tilde_qCM.E()+((n4Tilde_q*q4)+(n4Tilde_q*PBeam))*n4_qCM.E(),
               sqrt(W2), sqrt(q4.Mass2()+PBeam.Mass2()+2.0*(q4*PBeam)));
        printf("   X4_q*X4_q = %10.6f = %10.6f = X4_qCM*X4_qCM, Y4_q*q4 = %13.6e, n4_q*n4Tilde_q = %10.6f \n",
               X4_q*X4_q, X4_qCM*X4_qCM, Y4_q*q4, n4_q*n4Tilde_q);
    }
     */
    return iEvt;
}


#endif /* Deep_event_cxx */
