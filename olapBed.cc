
#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <limits>
#include <map>
#include "trees/interval_tree.h"
#include <getopt.h>
#include <argp.h>

using namespace std;

/*
NOTE1: not sure how to remember precision of score.  printf %f prints 8 by default, I can make
      I can make it 20, but then it adds 'random' numbers to fill out any number that didn't
      have 20.  storing as string would remember it, but then i can't add them... What if a file
      has numbers of varying precision, you have to remember on a 1-by-1 basis.  For now, just output
      the standard 8. BUT have to remember.
*/ 

class myIntv : public Interval {
public:
  myIntv(const unsigned long low,const unsigned long high, const double scr, const string id) :_low(low), _high(high),_scr(scr),_id(id) { }
  
  unsigned long GetLowPoint() const { return _low;}
  unsigned long GetHighPoint() const { return _high;}
  double GetScr() const { return _scr;}
  string GetId() const { return _id;} // not sure this is efficient... i want to keep this small, maybe there's another way
  IntervalTreeNode * GetNode() { return _node;}
  void SetNode(IntervalTreeNode * node) {_node = node;}
protected:
  unsigned long _low;
  unsigned long _high;
  double _scr;
  string _id;
  IntervalTreeNode * _node;
};

bool hasScore = 0;
bool dmyScore = 0;
string  outdlm("\t");

struct bed_s {
  string chr;
  unsigned long bgn;
  unsigned long end;
  string id;
  double scr;
  string strand;
  string otro;
};

typedef map<string, IntervalTree *> chromtree;

void print_bed(bed_s *b) {
  printf("%s\t%ld\t%ld\t%s",b->chr.c_str(),b->bgn-1,b->end, b->id.c_str());
  if (b->scr > -1*numeric_limits<double>::max()) printf("\t%f", b->scr); else if (dmyScore) printf("\t.");
  if (!b->strand.empty()) printf("\t%s", b->strand.c_str());
  if (!b->otro.empty()) printf("\t%s", b->otro.c_str());
}


chromtree fn2IntervalTree(string fn) {

  chromtree ctree;

  // ifstream fhb(fn.c_str()); if (!fhb.is_open()) {cerr<<"Unable to open file: "<<fn<<endl ;abort();}
  ifstream fbin;
  if (fn != "-") {
    fbin.open(fn.c_str());
    if (!fbin.is_open()) {cerr<<"Unable to open file: "<<fn<<endl ;abort();}
  }
  istream &fhb  = (fn=="-") ? std::cin : fbin;

  bed_s *cbed= new bed_s();

  //NOTE: need to error check that if parsing 
  while(fhb>>cbed->chr>>cbed->bgn>>cbed->end) {
    cbed->bgn++;
    // optional fields   
    cbed->id.clear();                         if (fhb.peek()!='\n') fhb>>cbed->id;
    
    // score could be .
    cbed->scr =-1*numeric_limits<double>::max(); 
    if (fhb.peek() !='\n'){
      fhb.get();
      if (fhb.peek() != '.') fhb>>cbed->scr; 
    }

    fhb.ignore(numeric_limits<streamsize>::max(),'\n');

       /* 
       printf("Load db.bed: %s %ld %ld [%s]",cbed->chr.c_str(),cbed->bgn-1,cbed->end, cbed->id.c_str());
       if (cbed->scr > numeric_limits<double>::min()) printf("  opt:%f %s %s", cbed->scr,cbed->strand.c_str(),cbed->otro.c_str());
       printf("\n");
       print_bed(cbed);printf("\n");
*/
       


    if (ctree.find(cbed->chr) == ctree.end()) ctree[cbed->chr] = new IntervalTree();
    ctree[cbed->chr]->Insert(new myIntv(cbed->bgn,cbed->end,cbed->scr,cbed->id));
  }

  if (!fhb.eof() ) cerr<<fn<<" err:"<<fhb.rdstate()<<" !!! db not all parsed!\n";
  //fhb.close();

  return ctree;
}


int main ( int argc, char *argv[] ) {

  /* extensions
     1) can do multiple qrys
     2) % olap                --done
     3) output each individual bedB or bedB scores --done
     4) output intersection
     5) 
     */

  int c;
  bool ocount   = 0;
  bool oscore   = 0;
  bool allscore = 0;
  bool allid    = 0;
  bool allolap   = 0;
  bool allline   =0;
  float db_olap  =0.0;
  float qry_olap = 0.0;
  while ( (c = getopt(argc, argv, "acsiobd:q:m:")) != -1)   //c: count olaps, s: score overlaps, a: print scr of overlaps
    switch (c) {
      case 'c':  ocount = 1; break;
      case 's':  oscore = 1; break;
      case 'a':  allscore = 1; break;
      case 'i':  allid = 1; break;
      case 'o':  allolap = 1;break;
      case 'b':  allline = 1;break;
      case 'd':  db_olap = atof(optarg);break;
      case 'q':  qry_olap = atof(optarg);break;
      case 'm':  outdlm.assign(optarg); break;
      case '?':
                 if (optopt == 'c')         fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                 else if (isprint (optopt)) fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                 else                       fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
                 return 1;
      default:
                 abort ();
    }
  if (optind+2 != argc ) {
    cout<<"usage: "<< argv[0]<<
      "[options] bedA bedB [bedA is db, bedB is qry.]\n"<<
      "Note: by default only those lines in bedB overlapping a line in bedA are reported\n"<<
      "Note: Right now, if bed score isn't numeric or ., parsing of the line stops there and screws things up -need to change maybe .\n"<<
      "Options\n"<<
      "-c  append the number of db lines olapping each qry line [including 0]\n"<<
      "-s  append the sum of db line scores to each qry line [including 0]\n"<<
      "-a  append all of the db line scores to each qry line [comma separated inluding none]\n"<< 
      "-i  append all of the db line ids to each qry line [comma separated inluding none]\n"<< 
      "-o  append all of the overlap between db and qry line [comma separated inluding none]\n"<< 
      "-b  append all of the db lines to each qry line [comma separated]\n"<< 
      "-m  delimiter for appended info [dflt: <tab>]\n"<<
      "-d  amount of olap of db req'd\n"<<
      "-q  amount of olap of q req'd\n"<<
      // just the intersection
      endl; 
    abort();
  } 
  string abed(argv[argc-2]); 
  string bbed(argv[argc-1]); 
  //printf("cnt:%d scr:%d abed:%s bbed:%s optind:%d %d\n",ocount,oscore,abed.c_str(),bbed.c_str(), optind,argc);
  if (abed == "-" && bbed=="-") {cout<<"both a and b bed can't be -\n";abort();}

  chromtree ctree = fn2IntervalTree(abed);

// ifstream fh(bbed.c_str()); if (!fh.is_open()) {cerr<<"Unable to open file: "<<bbed<<endl ;abort();}
  ifstream fbin;
  if (bbed != "-") {
    fbin.open(bbed.c_str());
    if (!fbin.is_open()) {cerr<<"Unable to open file: "<<bbed<<endl ;abort();}
  }
  istream &fh  = (bbed=="-") ? std::cin : fbin;



 fh.exceptions ( ifstream::badbit| ifstream::badbit );


  bed_s *cbed = new bed_s();
  chromtree::iterator tree;



  bool print_all = (ocount||oscore||allscore||allid||allolap||allline);
  while(fh >>cbed->chr >>cbed->bgn >>cbed->end) {

    cbed->bgn++;
    // optional fields   
    cbed->id.clear();                         if (fh.peek()!='\n') fh>>cbed->id;
    // score could be .
    cbed->scr =-1*numeric_limits<double>::max(); 
    if (fh.peek() != '\n'){
      hasScore=1; 
      fh.get();
      if (fh.peek() == '.'){ dmyScore=1; fh.get();}
      else                 { fh>>cbed->scr; }

   }
    cbed->strand.clear();                     if (fh.peek()!='\n') fh>>cbed->strand;
    cbed->otro.clear();                       if (fh.peek()!='\n') std::getline(fh, cbed->otro);

//chr1    10619   10621   .       1:1     :
//chr1    10788   10790   .       :       0:1


//  cout<<cbed->chr<<" "<<cbed->bgn<<" "<<cbed->end<<" "<<cbed->id<<" "<<cbed->scr<<" "<<cbed->strand<<" "<<cbed->otro<<endl;

//    std::cout.setf(std::ios::unitbuf); 
    if ( (tree = ctree.find(cbed->chr)) != ctree.end()) {
      // (Note: these aren't returned in order)
      TemplateStack<void *> *r = tree->second->Enumerate(cbed->bgn,cbed->end);
      
      if (qry_olap > 0.0 || db_olap > 0.0) {
        TemplateStack<void *> *tmp = new TemplateStack<void *>(r->Size());
        //cout<<qry_olap<<" "<<db_olap<<endl;        
        for (int j=0; j< r->Size(); j++)  {
          unsigned long ab = ((myIntv*)(*r)[j])->GetLowPoint();
          unsigned long ae = ((myIntv*)(*r)[j])->GetHighPoint();
          unsigned long olap = ((cbed->end > ae)? ae : cbed->end) 
            - ((cbed->bgn > ab)? cbed->bgn : ab) +1;
          if ( (db_olap == 0.0 || olap >= db_olap*(ae - ab)) && 
              (qry_olap == 0.0 || olap >= qry_olap*(cbed->end - cbed->bgn)) ) {
            //cout<<j<<" "<<olap<<">= d"<<(ae-ab)*db_olap<<" q"<<(cbed->end-cbed->bgn)*qry_olap<<" --- "<<ae<<"-"<<ab<<" "<<cbed->end<<"-"<<cbed->bgn<<endl;
            tmp->Push( (myIntv*)(*r)[j] );
          } 
        }
        r->Clear();
        r=tmp;
        //cout<<r->Size()<<" "<<tmp->Size();//return(0);
      }

      // print bedB line
      if (print_all || r->Size() > 0) print_bed(cbed); 

      // append results
      if (ocount) { printf("%s%d",outdlm.c_str(),r->Size()); } 
      if (oscore) {
        double score = 0;
        for (int j=0; j< r->Size(); j++) 
          if (((myIntv*)(*r)[j])->GetScr() > -1*numeric_limits<double>::max())
            score +=  ((myIntv*)(*r)[j])->GetScr();
        printf("%s%f",outdlm.c_str(),score);
      } 
      if (allscore) {
        cout<<outdlm; 
        int cnt = 0;
        for (int j=0; j< r->Size(); j++) 
          if (((myIntv*)(*r)[j])->GetScr() > -1*numeric_limits<double>::max()) {
            if (cnt++) printf(","); 
            printf("%f",((myIntv*)(*r)[j])->GetScr() );
          }
      }
      if (allolap) {
        cout<<outdlm; 
        int cnt = 0;
        for (int j=0; j< r->Size(); j++)  {
          if (cnt++) printf(","); 
          unsigned long ab = ((myIntv*)(*r)[j])->GetLowPoint();
          unsigned long ae = ((myIntv*)(*r)[j])->GetHighPoint();
          unsigned long olap = ((cbed->end > ae)? ae : cbed->end) 
                   - ((cbed->bgn > ab)? cbed->bgn : ab) +1;
          printf("%ld",olap);
        }
      }
      if (allid) {
        cout<<outdlm; 
        int cnt = 0;
        for ( long j=0; j< r->Size(); j++)  {
          if ( !((myIntv*)(*r)[j])->GetId().empty() ) {
            if (cnt++) printf(","); 
            printf("%s",((myIntv*)(*r)[j])->GetId().c_str() );
          }
        }
      }
      if (allline) {
        cout<<outdlm; 
        int cnt = 0;
          bed_s *rbed = new bed_s();
          for (int j=0; j< r->Size(); j++)  {
            if (cnt++) printf(","); 
            rbed->chr = cbed->chr;
            rbed->bgn = ((myIntv*)(*r)[j])->GetLowPoint();
            rbed->end = ((myIntv*)(*r)[j])->GetHighPoint();
            string id = ((myIntv*)(*r)[j])->GetId();
            if (!id.empty()) rbed->id = id;
            double scr = ((myIntv*)(*r)[j])->GetScr();
            if (scr > -1*numeric_limits<double>::max()) rbed->scr = scr;
            print_bed(rbed);
          }
       delete rbed;
      }
      if (print_all || r->Size() > 0) printf("\n");
      delete r; 
    //} else if (print_all) { print_bed(cbed); ocount ? printf("%s0\n", outdlm.c_str()) : printf("\n"); }
    } else if (print_all) { print_bed(cbed); ocount ? cout<<outdlm<<"0"<<endl : cout<<endl; }

  }

  if (!fh.eof()) cerr<<fh<<" err:"<<fh.rdstate()<<" !!! qry not all parsed! "<<endl;

//  fh.close(); 
}
