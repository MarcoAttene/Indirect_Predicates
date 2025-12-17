DEL indirect_predicates.h
..\x64\Release\jpck.exe > indirect_predicates.h
FOR %%f in (direct\*.txt) DO ..\x64\Release\jpck.exe %%f >> indirect_predicates.h
FOR %%f in (indirect\*.txt) DO ..\x64\Release\jpck.exe %%f >> indirect_predicates.h
copy /Y indirect_predicates.h ..\..\include\