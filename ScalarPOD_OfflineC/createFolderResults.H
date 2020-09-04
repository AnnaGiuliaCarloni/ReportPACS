fileName currentFolder= args.caseName ();

if (isMultiDir)
{
    if(isDir(args.rootPath()/"POD_"+fieldName))
    {
        Info<<"\n'POD_"<<fieldName<<"' folder already exists" << endl;
        Info<<"\n'deleting it" << endl;
        rmDir(args.rootPath()/"POD_"+fieldName);
    }

    Info<<"\n creating ''POD_"<<fieldName<<"'' folder \n" << endl;
    mkDir(args.rootPath()/"POD_"+fieldName);
    mkDir(args.rootPath()/"POD_"+fieldName/"0");
    cp("system", args.rootPath()/"POD_"+fieldName/"system");
    cp("constant", args.rootPath()/"POD_"+fieldName/"constant");


    chDir(args.rootPath()/"POD_"+fieldName);
    currentFolder="POD_"+fieldName;
}


Foam::Time runTime(Foam::Time::controlDictName, args.rootPath(), currentFolder);
runTime.setTime(0,0);
