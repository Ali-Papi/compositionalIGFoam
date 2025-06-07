{
    fvScalarMatrix xCO2Eqn
        (
            eps*fvm::ddt(Sb,xCO2) + fvm::div(phib,xCO2)
           // ==
            // event source terms
           // - sourceTerm
        );
    xCO2Eqn.solve();
    
    xH2O = 1-xCO2;
    
    // Flash:
    //NL = NH2O/(1-PH) - NCO2/(PP-1);
    //NG = NH2O + NCO2 - NL;
    //NH2Ol = xH2O*NL;
    //NCO2l = xCO2*NL;
    //NH2Og = yH2O*NG;
    //NCO2g = yCO2*NG;
    Info << "xCO2: " << xCO2 << endl;
    Info << "xCO2: " << " Min(xCO2) = " << gMin(xCO2) << " Max(xCO2) = " << gMax(xCO2) << endl;

}
