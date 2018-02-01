REM Write a  message for this commit.
set /p commit_message=Commit message:

git add -A
git commit -m "%commit_message%"
git push
