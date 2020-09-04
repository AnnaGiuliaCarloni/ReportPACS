/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    ScalarPOD_T

Description

\*---------------------------------------------------------------------------*/
#include <iostream>
#include <fstream>
#include <fstream>
#include <iomanip>

#include "fvCFD.H"
#include "argList.H"
#include "timeSelector.H"
#include "volFields.H"
#include "MOR.H"
#include "IOmanip.H"
#include "IOstream.H"
#include "POD_EigenBase.H"

// * * * * * * * * * * * * * * * *LOCAL INCLUDES* * * * * * * * * * * * * * * //

#include "getPODOfflineParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    auto start_s=std::clock();

#include "addDictOption.H"
    timeSelector::addOptions(false, true);
#include "setRootCase.H"

    //get POD parameters

    PODparameters POD_parameters = getPODparameters(args);

    List<fileName> folders_list = POD_parameters.folders_list ;

    word fieldName ("T");

    label maxBasis = POD_parameters.maxBasis;

    scalar accuracy = POD_parameters.accuracy;

    // Initialise all Lists necessary to store the snapshots and PODBasis
    // I used a list because I have to store data of type volScalarField

    PtrList<volScalarField> scalarSnapshotsList;
    PtrList<volScalarField> scalarPODOrthoNormalBasis;

    // get ScalarSnapshots

    for (label folderI=0; folderI < folders_list.size(); ++ folderI)
    {
        chDir(args.rootPath()/folders_list[folderI]);

        Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), folders_list[folderI]);

#include "CreateMesh.H"

        Info << "reading snapshots in "<< folders_list[folderI]<<"\n"<< endl;

        // Get times list
        instantList timeDirs = timeSelector::select0(runTime, args);

        forAll(timeDirs, timeI)
        {
            runTime.setTime(timeDirs[timeI], timeI);

            Info<< "Time = " << runTime.timeName() << endl;

            if (mesh.readUpdate() != polyMesh::UNCHANGED)
            {

                FatalErrorInFunction
                        << "polyMesh has changed. POD can be performed only on unchanged polyMesh"
                        << abort(FatalError);
            }

            scalarSnapshotsList.append
            (
                new volScalarField

                (
                    IOobject
                    (
                        fieldName,
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ
                    ),
                    mesh
                )
            );
        }

    }

    // create new folder in which all POD Basis will be stored

#include"createFolderResults.H"

#include "CreateMesh.H"

    // performPOD

    POD_EigenBase eigenBase( scalarSnapshotsList );

    label Basisize = 0;

    const scalarField& EigenValues = eigenBase.eigenValues();
    const scalarField& cumEigenValues = eigenBase.cumulativeEigenValues();

    forAll (cumEigenValues, i)
    {
        Basisize++;

        if ((cumEigenValues[i] > accuracy) | (Basisize==maxBasis))
        {
            break;
        }
    }

    Info<< "Cumulative eigen-values: "
        << setprecision(14) << cumEigenValues << nl
        << "Base size: " << Basisize << endl;

    /// Establish orthonormal base

    scalarPODOrthoNormalBasis.setSize(Basisize);

    for (label baseI = 0; baseI < Basisize; baseI++)
    {
        const scalarField& eigenVector = eigenBase.eigenVectors()[baseI];

        GeometricField<scalar, fvPatchField, volMesh> * onBasePtr
        (
            new GeometricField<scalar, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldName + "POD" + name(baseI),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ
                ),
                mesh,
                dimensioned<scalar>
                (
                    "zero",
                    scalarSnapshotsList[baseI].dimensions(),
                    pTraits<scalar>::zero
                )
            )
        );

        GeometricField<scalar, fvPatchField, volMesh>& onBase = *onBasePtr;

        forAll (eigenVector, eigenI)
        {
            onBase += eigenVector[eigenI]*scalarSnapshotsList[eigenI];
        }

        /// Re-normalise ortho-normal vector

        scalar SqrtEigenValue = Foam::sqrt(EigenValues[baseI]);
        if (SqrtEigenValue > SMALL)
        {
            onBase /= SqrtEigenValue;
        }

        scalarPODOrthoNormalBasis.set(baseI, onBasePtr);
    }

    // Calculate interpolation coefficients(io non ho bisogno di usare l'interpolazione) 

    scalarRectangularMatrix coeffs(scalarSnapshotsList.size(), scalarPODOrthoNormalBasis.size());

    forAll (scalarSnapshotsList, snapshotI)
    {
        forAll (scalarPODOrthoNormalBasis, baseI)
        {
            coeffs[snapshotI][baseI] =
                MOR::projection
                (
                    scalarSnapshotsList[snapshotI],
                    scalarPODOrthoNormalBasis[baseI]
                );
        }
    }

    // Check all snapshots

    forAll (scalarSnapshotsList, snapI)
    {

        volScalarField reconstruct
        (
            IOobject
            (
                fieldName + "PODreconstruct",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                scalarSnapshotsList[snapI].dimensions(),
                pTraits<scalar>::zero
            )
        );

        for (label baseI = 0; baseI < scalarPODOrthoNormalBasis.size(); baseI++)
        {
            reconstruct +=
                coeffs[snapI][baseI]*scalarPODOrthoNormalBasis[baseI];
        }


        scalar sumFieldError =
            MOR::L2norm
            (
                reconstruct - scalarSnapshotsList[snapI]
            );

        scalar measure =
            MOR::L2norm(scalarSnapshotsList[snapI]) + SMALL;

        Info<< "Field error: absolute = " << sumFieldError << " relative = "
            << sumFieldError/measure << " measure = " << measure
            << endl;

    }

    // Write out all fields

    Info << "Writing POD scalar base for volScalarField "
         << fieldName << " in folder '"
         << runTime.timeName()<<"'"  << endl;

    for (label i = 0; i < scalarPODOrthoNormalBasis.size(); i++)
    {
        scalarPODOrthoNormalBasis(i)->write();
    }

    

//Print out all the matrices

scalarRectangularMatrix L (scalarPODOrthoNormalBasis.size(),scalarPODOrthoNormalBasis.size(), Foam::zero());

for (label ii=0; ii< scalarPODOrthoNormalBasis.size(); ii++)
{
	for (label jj=0; jj< scalarPODOrthoNormalBasis.size(); jj++)
        {
		L[ii][jj]=fvc::domainIntegrate(scalarPODOrthoNormalBasis[jj]*(fvc::laplacian(scalarPODOrthoNormalBasis[ii]))).value();
	}

}


std::ofstream f1 ("Lmatrix.txt");

if (f1.is_open())
{
	for (label ii=0; ii<scalarPODOrthoNormalBasis.size(); ii++)
	{
		for (label jj=0; jj<scalarPODOrthoNormalBasis.size(); jj++ )
		{
			f1 << std::setprecision(14) <<L(ii,jj)<< " ";
		}
		f1 << "\n";
	}
	f1.close();	
}


//Info << "Matrix L is " << L <<endl;



scalarRectangularMatrix C (scalarPODOrthoNormalBasis.size(),1, Foam::zero());

for (label ii=0; ii< scalarPODOrthoNormalBasis.size(); ii++)
{	        
		C[ii][0]=fvc::domainIntegrate(scalarPODOrthoNormalBasis[ii]).value();
}


std::ofstream f2 ("Cmatrix.txt");

if (f2.is_open())
{
	for (label ii=0; ii<scalarPODOrthoNormalBasis.size(); ii++)
	{
		f2 << std::setprecision(14) <<C(ii,0)<< " ";
		f2 << "\n";
	}
	f2.close();	
}

//Info << "Matrix C is " << C <<endl; 


scalarRectangularMatrix D (scalarPODOrthoNormalBasis.size(),1, Foam::zero());

for (label ii=0; ii< scalarPODOrthoNormalBasis.size(); ii++)
{	        
                scalar SurfaceIntegral=0;
		forAll (mesh.boundary()[mesh.boundary().findPatchID("atmosphere")].magSf(), SfI)
		{
			SurfaceIntegral+=((mesh.boundary()[mesh.boundary().findPatchID("atmosphere")].magSf()[SfI])*(scalarPODOrthoNormalBasis[ii].boundaryField()[ mesh.boundary().findPatchID("atmosphere")][SfI]));
		}
		D[ii][0]=SurfaceIntegral;

}


std::ofstream f3 ("Dmatrix.txt");

if (f3.is_open())
{
	for (label ii=0; ii<scalarPODOrthoNormalBasis.size(); ii++)
	{
		f3 << std::setprecision(14) <<D(ii,0)<< " ";
		f3 << "\n";
	}
	f3.close();	
}

//Info << "Matrix D is " << D <<endl;

scalarRectangularMatrix B (scalarPODOrthoNormalBasis.size(),scalarPODOrthoNormalBasis.size(), Foam::zero());

for (label ii=0; ii< scalarPODOrthoNormalBasis.size(); ii++)
{
	for (label jj=0; jj< scalarPODOrthoNormalBasis.size(); jj++)
        {
                scalar SurfaceIntegral=0;
		forAll (mesh.boundary()[mesh.boundary().findPatchID("atmosphere")].magSf(), SfI)
		{
			SurfaceIntegral+=((mesh.boundary()[mesh.boundary().findPatchID("atmosphere")].magSf()[SfI])*(scalarPODOrthoNormalBasis[jj].boundaryField()[ mesh.boundary().findPatchID("atmosphere")][SfI])*(scalarPODOrthoNormalBasis[ii].boundaryField()[ mesh.boundary().findPatchID("atmosphere")][SfI]));
		}
		B[ii][jj]=SurfaceIntegral;  //faccio integrale di Diff a partee poi moltiplico per ciascuna componente
                                            //opuure lo metto tutto nello stesso integrale?
	}

}


std::ofstream f4 ("Bmatrix.txt");

if (f4.is_open())
{
	for (label ii=0; ii<scalarPODOrthoNormalBasis.size(); ii++)
	{
		for (label jj=0; jj<scalarPODOrthoNormalBasis.size(); jj++ )
		{
			f4 << std::setprecision(14) <<B(ii,jj)<< " ";
		}
		f4 << "\n";
	}
	f4.close();	
}

//Aggiunta per calcolo condizioni iniziali coeffs (CIC)

scalarRectangularMatrix CIC (scalarPODOrthoNormalBasis.size(), 1, Foam::zero());

forAll (scalarPODOrthoNormalBasis, basisIndex)
{
	CIC[basisIndex][0]=fvc::domainIntegrate(scalarPODOrthoNormalBasis[basisIndex]*scalarSnapshotsList[0]).value(); //solo con withZero
}

std::ofstream f8 ("CIC.txt");

if (f8.is_open())
{
	for (label ii=0; ii<scalarPODOrthoNormalBasis.size(); ii++)
	{
		
		
			f8 << std::setprecision(14) <<CIC(ii,0)<< " ";
		
		        f8 << "\n";
	}
	f8.close();	
}


#include "printAllFiles.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    auto stop_s=std::clock();

    Info<< nl << "ExecutionTime = " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
