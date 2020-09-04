
	    Dx=(D1_0*exp(-E1_0/T)+D2_0*sqrt(frate)*exp(-E2_0/T)+D3_0*frate) / Foam::pow(R.value(), 2);
            Dy=(D1_0*exp(-E1_0/T)+D2_0*sqrt(frate)*exp(-E2_0/T)+D3_0*frate)  / Foam::pow(R.value(), 2);
            Dz=(D1_0*exp(-E1_0/T)+D2_0*sqrt(frate)*exp(-E2_0/T)+D3_0*frate)  / Foam::pow(L.value(), 2);

volScalarField Diffxx = Dx;
volScalarField Diffxy = 0.0*Dx;
volScalarField Diffxz = 0.0*Dx;

volScalarField Diffyx = 0.0*Dy;
volScalarField Diffyy = Dy;
volScalarField Diffyz = 0.0*Dy;

volScalarField Diffzx = 0.0*Dz;
volScalarField Diffzy = 0.0*Dz;
volScalarField Diffzz = Dz;

forAll (Diff, cellI)
{
	Diff[cellI] = tensor(Diffxx[cellI],  Diffxy[cellI],  Diffxz[cellI],  Diffyx[cellI],  Diffyy[cellI],  Diffyz[cellI],  Diffzx[cellI],  Diffzy[cellI],  Diffzz[cellI]);
}

Diff.correctBoundaryConditions();

        C.correctBoundaryConditions();
        {
           

            fvScalarMatrix cEqn
            (
                Foam::pow(L.value(),2) * fvm::ddt(C) - Foam::pow(L, 2)*fvm::laplacian(Diff, C) == frate*gasSource * Foam::pow(L.value(),2)
            );

            cEqn.solve();
		}


