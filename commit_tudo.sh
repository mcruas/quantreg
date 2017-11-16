git add .
DATE=`date '+%d-%m-%Y'`
msg_commit='Commit do dia '$DATE
git commit -m "$msg_commit"
git push
