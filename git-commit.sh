#!/bin/bash
#git status
echo 'What changes have you made?'
read changes && export changes
echo 'What files would you like to commit? (Type -A for all of them)'
read files

git add $files
git status
git commit -m "$changes"
git pull
git push

echo "${GREEN}Thanks for improving our code!${NC}"
