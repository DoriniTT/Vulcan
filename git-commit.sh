#!/bin/bash
#git status
echo 'What changes have you made?'
read changes && export changes
echo 'What files would you like to commit? (Type -A for all of them)'
read files

git add $files
git status
git commit -m "$changes"
git pull origin master
git push origin master

GREEN='\033[0;32m'
NC='\033[0m'

echo "Would you like to send an email to the collaborators? (y/n)"
read answer && [ $answer = "y" ] && python send_email.py \
            && echo "Please set your email configurations in the 'send_email.py' script."

echo "${GREEN}Thanks for improving our code!${NC}"
