Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading diffusivity DT\n" << endl;

dimensionedScalar DT
(
    transportProperties.lookup("DT")
);

Info<< "Reading fission Rate\n" << endl;

dimensionedScalar frate
(
    transportProperties.lookup("frate")
);

Info<< "Reading field C\n" << endl;

volScalarField C
(
    IOobject
    (
        "C",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


dimensionedScalar D1_0
(
	"D1_0",
	transportProperties.lookup("D1_0")
);

dimensionedScalar E1_0
(
	"E1_0",
	transportProperties.lookup("E1_0")
);

dimensionedScalar D2_0
(
	"D2_0",
	transportProperties.lookup("D2_0")
);

dimensionedScalar E2_0
(
	"E2_0",
	transportProperties.lookup("E2_0")
);
dimensionedScalar D3_0
(
	"D3_0",
	transportProperties.lookup("D3_0")
);
dimensionedScalar L
(
	"L",
	transportProperties.lookup("L")
);
dimensionedScalar R
(
	"R",
	transportProperties.lookup("R")
);

Info<< "Reading diffusivity Diff\n" << endl;


volScalarField Dx
(
	IOobject
	(
		"Dx",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), pTraits<scalar>::zero)
);
volScalarField Dy
(
	IOobject
	(
		"Dy",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), pTraits<scalar>::zero)
);
volScalarField Dz
(
	IOobject
	(
		"Dz",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedScalar("zero",dimensionSet(0,0,-1,0,0,0,0), pTraits<scalar>::zero)
);

//// VECCHIA DEFINIZIONE DI TENSORE SIMMETRICO COMMENTATA
/*
volSymmTensorField Diff
(
	IOobject
	(
		"Diff",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
    mesh,
	dimensionedSymmTensor
    ( 
	   "Diff",
        dimensionSet(0, 2, -1, 0, 0, 0, 0),
	pTraits<symmTensor>::zero    //  
    ),
    zeroGradientFvPatchTensorField::typeName//SIMO
    
);*/

///// added by Eric
volTensorField Diff
(
    IOobject
    (
        "Diff",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
dimensionedTensor
    ( 
	   "Diff",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        pTraits<tensor>::zero
    )
);

volScalarField heatSource
	(
		IOobject
		(
			"heatSource",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("", dimensionSet(0,0,-1,1,0,0,0),7.27e-18),
		zeroGradientFvPatchScalarField::typeName
	);

volScalarField gasSource
	(
		IOobject
		(
			"gasSource",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("", dimensionSet(0,0,-1,0,0,0,0),5e-25),
		zeroGradientFvPatchScalarField::typeName
	);

//scalarSquareMatrix temp_coeffsMatrix(IFstream(Documenti/filedname+_coeffs_out.txt)());






