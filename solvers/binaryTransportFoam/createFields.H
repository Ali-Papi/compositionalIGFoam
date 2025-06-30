///////////////////////////////////////////////////////////////
////////////////////// TIMESTEP CONTROL ///////////////////////
///////////////////////////////////////////////////////////////

word CFL =  runTime.controlDict().lookupOrDefault<word>("CFL", "Coats");
if (adjustTimeStep) adjustTimeStep=true; // to remove warnings at compilation
int CFLint = -1;
if (CFL == "Coats") CFLint = 0;
else if (CFL == "Todd") CFLint = 1;
else if (CFL == "Courant") CFLint = 2; 
else
{
    FatalErrorIn
        (
            "in createFields.H"
        )
        << "CFL condition unrecongnized : Coats, Todd and Courant available" 
            << exit(FatalError);
}

//////////////////////////////////////////////////////////////////
////////////////////// PRESSURE SATURATION ///////////////////////
//////////////////////////////////////////////////////////////////

Info << "Reading pressure field p" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading saturation field Sb" << endl;
volScalarField Sb
(
    IOobject
    (
        "Sb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading CO2 field xCO2" << endl;
volScalarField xCO2
(
    IOobject
    (
        "xCO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info << nl << "Reading CO2 field yCO2" << endl;
volScalarField yCO2
(
    IOobject
    (
        "yCO2",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info << nl << "Reading H2O field xH2O" << endl;
volScalarField xH2O
(
    IOobject
    (
        "xH2O",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
   //dimensionedScalar("",dimless,1)
);
//xH2O=1-xCO2;

Info << nl << "Reading H2O field yH2O" << endl;

volScalarField yH2O
(
    IOobject
    (
        "yH2O",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //dimensionedScalar("",dimless,1)
);
//yH2O=1-yCO2;


scalar dSmax(runTime.controlDict().lookupOrDefault<scalar>("dSmax",0.));

//////////////////////////////////////////////////////////////////
////////////////////// THERMOPHYSICAL PROPERTIES /////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading thermophysicalProperties" << endl;

IOdictionary thermophysicalProperties
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dictionary subDict0 = thermophysicalProperties.subDict("solutionSpecies");
dictionary subDict1 = subDict0.subDict("Species");
//dimensionedScalar D1 = subDict1.lookup("D1");
//dimensionedScalar D2 = subDict1.lookup("D2"); //both phases should have it, if any.
dimensionedScalar H = subDict1.lookup("H");
dimensionedScalar Mw = subDict1.lookup("Mw")*1e-3;
dimensionedScalar Mg = subDict1.lookup("Mg")*1e-3;
dimensionedScalar Psat = subDict1.lookup("Psat");

volScalarField PH = p/H;
volScalarField PP = p/Psat;

//////////////////////////////////////////////////////////////////
////////////////////// TRANSPORT PROPERTIES //////////////////////
//////////////////////////////////////////////////////////////////

Info << nl << "Reading transportProperties" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

//- list that receives event files of event-based boundary conditions
List<patchEventFile*> patchEventList;
eventInfiltration::setEventFileRegistry(&patchEventList, Sb.name());

dimensionedScalar Sbmin(transportProperties.lookupOrDefault("Sbmin",dimensionedScalar("Sbmin",dimless,0)));
//Info << "Sbmin from extrction method 1 is " << Sbmin << endl;
//dimensionedScalar Sbmin2 = transportProperties.subDict("BrooksAndCoreyCoeffs").lookupOrDefault<dimensionedScalar>("Sbmin", dimensionedScalar("Sbmin", dimless, 0));
scalar Sbmin1 = transportProperties.subDict("BrooksAndCoreyCoeffs").lookupOrDefault<scalar>("Sbmin", 0);
//Info << "Sbmin from extrction method 2 is " << Sbmin1 << endl;
//scalar Sbmin = transportProperties.subDict("BrooksAndCoreyCoeffs").lookupOrDefault<scalar>("Sbmin", 0.0);
dimensionedScalar SbminDim("Sbmin", dimless, Sbmin1);
scalar Sbmax1 = transportProperties.subDict("BrooksAndCoreyCoeffs").lookupOrDefault<scalar>("Sbmax", 0);
dimensionedScalar SbmaxDim("Sbmin", dimless, Sbmax1);

/////////////////////////////////////////////////////////////////////////////
/////////////////////////// PHASE MODEL CREATION ////////////////////////////
/////////////////////////////////////////////////////////////////////////////

autoPtr<incompressiblePhase> phasea = incompressiblePhase::New(mesh,transportProperties,"a");
volVectorField& Ua = phasea->U();
surfaceScalarField& phia = phasea->phi();
const dimensionedScalar& rhoa = phasea->rho();
const dimensionedScalar& mua = phasea->mu();

autoPtr<incompressiblePhase> phaseb = incompressiblePhase::New(mesh,transportProperties,"b");
volVectorField& Ub = phaseb->U();
surfaceScalarField& phib = phaseb->phi();
const dimensionedScalar& rhob = phaseb->rho();
const dimensionedScalar& mub = phaseb->mu();

dimensionedScalar Mmu(mub/mua);//viscosity ratio

/////////////////////////////////////////////////////////////////////////////
////////////////////////// POROUS MEDIA PROPERTIES //////////////////////////
/////////////////////////////////////////////////////////////////////////////

// Relative permeability / Capillary pressure models 
autoPtr<twophasePorousMediumModel> porousModel = twophasePorousMediumModel::New("Sb", mesh, transportProperties, phaseb);
autoPtr<relativePermeabilityModel>& krModel = porousModel->krModel();
autoPtr<capillarityModel>& pcModel = porousModel->pcModel();
scalar activateCapillarity(transportProperties.lookupOrDefault<scalar>("activateCapillarity",0.));

// Porosity
porousModel->check_eps();
const volScalarField& eps = porousModel->eps();

// Intrinsic permeability
porousModel->check_K();
const volScalarField& K = porousModel->K();

// permeability interpolation
surfaceScalarField Kf(fvc::interpolate(K,"K"));

/////////////////////////////////////////////////////////////////////////////
////////////////////////// VELOCITY - FLUXES ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

Info << nl << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    Ua + Ub
);

#include "createPhi.H"
surfaceScalarField phiP = phi;

///////////////////////////////////////////////////////////////////
////////////////////////// FORCING TERMS //////////////////////////
///////////////////////////////////////////////////////////////////

volScalarField sourceTerm
(
    IOobject
    (
        "sourceTerm",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimless/dimTime,0)
);

double Sb_error = 1;
volScalarField SbLastItr = Sb;
const volScalarField Sb0 = Sb;
volScalarField SbLastItr1 = Sb;
volScalarField NH2O = eps*(xH2O*rhob*Sb/Mw+yH2O*rhoa*(1-Sb)/Mg);
volScalarField NCO2 = eps*(xCO2*rhob*Sb/Mw+yCO2*rhoa*(1-Sb)/Mg);
volScalarField NT = eps*(rhob*Sb/Mw+rhoa*(1-Sb)/Mg);
volScalarField NL = eps*rhob*Sb/(Mw);
volScalarField NG = NT-NL;
volScalarField NLmax = eps*rhob*SbmaxDim/Mw; //Sb here means Sbmax in Pc equation
volScalarField NLmin = eps*rhob*SbminDim/Mw; //Sb here means Sbmin in Pc equation, for example 0.9999
volScalarField NH2Ol = eps*xH2O*rhob*Sb/(Mw);
volScalarField NH2Og = eps*yH2O*rhoa*(1-Sb)/(Mg);
volScalarField NCO2l = eps*xCO2*rhob*Sb/(Mw);
volScalarField NCO2g = eps*yCO2*rhoa*(1-Sb)/(Mg);
volScalarField NCO2g_t = NCO2g;
volScalarField NH2Ol_t = NH2Ol;
volScalarField delNH2Ol = NH2Ol_t - NH2Ol;
volScalarField delNCO2g = NCO2g_t - NCO2g;
volScalarField Nnet_change = -delNH2Ol + delNCO2g;
volScalarField NLflash = NL;