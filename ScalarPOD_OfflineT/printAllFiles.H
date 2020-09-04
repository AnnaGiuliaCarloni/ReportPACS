// print files.txt

Info << "\nprinting files"<< endl;

if (isDir(args.rootPath()/currentFolder/fieldName+"_POD_files"))
{
    rmDir(args.rootPath()/currentFolder/fieldName+"_POD_files");
}

mkDir(args.rootPath()/currentFolder/fieldName+"_POD_files");
chDir(args.rootPath()/currentFolder/fieldName+"_POD_files");


ofstream EigenValuesfile (fieldName+"EigenValues.txt", std::ios::out | std::ios::app );
if (EigenValuesfile.is_open())
{
    for (label jj=0; jj< EigenValues.size(); jj++)
    {
        EigenValuesfile << std::setprecision(14)<< EigenValues[jj] ;
        EigenValuesfile << "\n";

    }
    EigenValuesfile.close();
}



