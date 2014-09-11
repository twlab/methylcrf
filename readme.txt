

# re-connect olapBed as a subtree
# Nonetheless, don't update olapBed here.  Clone the olapBed repo and
# make changes there -it's easier to keep a consistent history.
git remote add -f olapbed https://mscse@bitbucket.org/mscse/olapbed.git                                                                                   
git merge -s ours --no-commit --squash olapbed/master                                                                                                  
git pull -s subtree olapbed master

