//k=1/2(u'u'-v'v'-w'w')
//u' = u-u_mean
// volScalarField kRes = 0.5*(r(db.lookupObject<volSymmTensorField>("UPrime2Mean")); //Instantaneous resolved TKE
volVectorField UPrime = (u - db.lookupObject<volVectorField>("UMean")); //U' = u-u_mean
volScalarField kRes = 0.5*(UPrime & UPrime); //Instantaneous resolved TKE

volScalarField kSGS = turbulence->k();

volScalarField Ktot = kRes + kSGS;

LESResIndex = kRes/max(Ktot, kSmall);

Info << "ESTIMATED LESResIndex min, max, avg:" << min(LESResIndex).value()
    << ", " << max(LESResIndex).value()
    << ", " << LESResIndex.weightedAverage(mesh.V()).value() << endl;

Warning << "Averaging of velocity U is not active in the controlDict." << endl
    << "This is necessary to calculate resolved TKE" << endl
    << "Add field U to the fieldAverage utility in controlDict." << endl;


/////////////////////////////////////////////////////////////
void scalarField LESResIndex
{
    IObject
    {
    "LESResIndex",
    runTime.timeName(),
    mesh,
    IObject::READ_IF_PRESENT,
    IObject::AUTO_WRITE
    },
    mesh,
    dimensionedScalar
    {
    "LesResIndex",
    dimless,
    0.0
    }
};



const volScalarField kSmall
{
    IOobject
    {
    "kSmall",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
    },
    mesh,
    dimensionedScalar
    {
    "kSmall",
    dimVelocity*dimVelocity,
    SMALL
    }
};



//////////////CONTRIL DICT?????????????????????????????

functions
{
    FA1
    {
        functionObjectLibs ("libfieldFunctionObjects.so");
        fields
        {
            U
            {
                mean        on;
                prime2Mean  on;
                base        time;
            }
        };
        executeControl    timestep;
        writeControl      writeTime;
        type             fieldAverage;
        startTime        0;
    }
}

///////////////turbulence model///////////////

FoamFile
{
    version     2.0;
    format     ascii;
    class      dictionary;
    object     turbulenceProperties;
}
// ************************************************************** //

simulationType    LES;

LES
{
    LESModel    dynamicKEqn;
    turbulence    on;
    printCoeffs    on;
    delta    cubeRootVol;
    cubeRootVolCoeffs
    {
        deltaCoeff    1;
    }
}



// ************************************************************** //


boundaryField
{
    inlet
    {
        type            turbulentInlet;
        referenceField  uniform (1 0 0);
        fluctuationsScale (0.002 0.001 0.001);
        value           uniform (1 0 0);
    }
}

// ************************************************************** //

volScalarField Uwindowed
(
    IOobject
    (
        "Uwindowed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "Uwindowed", 
        dimVelocity,
        0.0
    )
);



IOobject dictIO
(
    "windowProperties",
    mesh.time().constant(),
    mesh,
    IOobject::MUST_READ
);

dictionary windowProperties;
windowProperties = IOdictionary(dictIO);

vector extents_min = windowProperties.lookupOrDefault<vector>("extents_min", vector(0,0,0));
vector extents_max = windowProperties.lookupOrDefault<vector>("extents_max", vector(100,100,100));

Info << "Extents min = " << extents_min << " and max = " << extents_max << "\n" << endl;

DynamicList<scalar> cellExtents;

forAll(mesh.C(), cellI)
{
    if (mesh.C()[cellI][0] > extents_min[0] && mesh.C()[cellI][0] < extents_max[0])
    {
        if (mesh.C()[cellI][1] > extents_min[1] && mesh.C()[cellI][1] < extents_max[1])
        {
            if (mesh.C()[cellI][2] > extents_min[2] && mesh.C()[cellI][2] < extents_max[2])
            {
                cellExtents.append(mag(mesh.C()[cellI]));
            }
        }
    }
}

Info<< "Total number of cells in the window: " << cellExtents.size() << "/" << mesh.C().size()<<endl;


scalar windowPeriod = windowProperties.lookupOrDefault<scalar>("window_period", 10);
scalar deltaT = runTime.deltaT().value();
scalar window = int(windowPeriod/deltaT);

Info << "Window Period = " << windowPeriod << "s" << nl << endl;
Info << "Number of entries (Window/deltaT) = " << window << nl << endl;



///////////////////////////////////////

List<FIFOStack>> velocityStacks(mesh.cells().size());
List<DynamicList<scalar>> velocityWindowed(mesh.cells().size());

/////////////////////////////////////////


forAll(cellExtents,label)
{
    int cellI = cellExtents[label];
    velocityStacks[cellI].push(mag(U[cellI]);


    if(cellI==cellExtents[0])
    {
    Info<< "velocityStacks after Push = " << velocityStacks[cellI] << nl << endl;
    }

    
    if (velocityStacks[cellI].size() == window)
    {
    scalar meanValue = 0;
    forAllConstIter(FIFOStack<scalar>, velocityStacks[cellI].iter)
    {
    meanValue += iter();
    }
    meanValue /= window;
    velocityStacks[cellI].pop();
    velocityWindow[cellI].append(meanValue);
    }
    Uwindowed[cellI] = velocityWindow[cellI].back();

    if(cellI==cellExtents[0])
    {
    Info<< "velocityStacks after Pop = " << velocityStacks[cellI] << nl << endl;
    Info<< "velocityWindow + meanValue " << velocityWindow[cellI] << nl << endl;
    Info<< "Uwindowed = " << Uwindowed[cellI] << nl << endl;
    }
    }
}