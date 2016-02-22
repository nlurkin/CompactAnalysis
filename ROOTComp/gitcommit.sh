#!/bin/sh

exit
currBranch=`git branch -v | grep "*" | cut -f 2 -d " "`

git checkout automatic
git pull origin automatic
git add CMakeLists.txt SelectionFunctions.* pi0d_Selection.cpp pid.cpp
git add CompactRead/*.h CompactRead/*.cpp CompactRead/CMakeLists.txt
git add Support/*.h Support/*.cpp Support/CMakeLists.txt
git commit -m "Automatic commit"
git push origin automatic

#autoCommit=`git branch -v | grep automatic | tr -s " " | cut -f 3 -d " "`
git reset $currBranch
git checkout $currBranch
#git cherry-pick --no-commit $autoCommit
