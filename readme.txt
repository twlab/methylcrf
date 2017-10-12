Updated on Dec/2016:
replace olapBed with mapBed from bedtools

——————————————————————

# NOTE for previous repo hook to tag every update
# Need to redo this functionality for general reproducibility
#  for users or think of another way
on my repo, periodically run:
git fetch --tags 
to get auto-updated tags

# re-connect olapBed as a subtree
# Nonetheless, don't update olapBed here.  Clone the olapBed repo and
# make changes there -it's easier to keep a consistent history.
git remote add -f olapbed https://mscse@bitbucket.org/mscse/olapbed.git                                                                                   
git merge -s ours --no-commit --squash olapbed/master                                                                                                  
git pull -s subtree olapbed master

