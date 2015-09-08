#!/bin/sh

currBranch=`git branch -v | grep "*" | cut -f 2 -d " "`

git checkout automatic
git pull origin automatic
git add CMakeLists.txt
git add batch/*.cpp
git add common/CMakeLists.txt common/*.h common/*.cpp
git add python/*.py
git add root/*.C
git commit -m "Automatic commit"
git push origin automatic

autoCommit=`git branch -v | grep automatic | tr -s " " | cut -f 3 -d " "`
git checkout $currBranch
git cherry-pick --no-commit $autoCommit
