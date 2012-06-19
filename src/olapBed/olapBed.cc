
#include <iostream>
#include <fstream>
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
  myIntv(const int low,const int high, const double scr, const string id) :_low(low), _high(high),_scr(scr),_id(id) { }
  
  int GetLowPoint() const { return _low;}
  int GetHighPoint() const { return _high;}
  double GetScr() const { return _scr;}
  string GetId() const { return _id;} // not sure this is efficient... i want to keep this small, maybe there's another way
  IntervalTreeNode * GetNode() { return _node;}
  void SetNode(IntervalTreeNode * node) {_node = node;}
protected:
  int _low;
  int _high;
  double _scr;
  string _id;
  IntervalTreeNode * _node;
};

bool hasScore = 0;
bool dmyScore = 0;
// bed struct
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

chromtree fn2IntervalTree(string fn) {

  chromtree ctree;

 ifstream fhb(fn.c_str()); if (!fhb.is_open()) {cerr<<"Unable to open file: "<<fn<<endl ;abort();}

  bed_s *cbed= new bed_s();

 //NOTE: need to error check that if parsing 
  string dummy;
  string dmyScore;
  while(fhb>>cbed->chr>>cbed->bgn>>cbed->end) {
     cbed->bgn++;
    // optional fields   
    cbed->id.clear();                         if (fhb.peek()!='\n') fhb>>cbed->id;
    // score could be .
    cbed->scr =-1*numeric_limits<double>::max(); 
    if (fhb.peek() !='\n'){hasScore=1; char d =fhb.get();
      if (fhb.peek() != '.'){ fhb.putback(d); fhb>>cbed->scr; }
      else                 { dmyScore=1;fhb.get();}
    }
   // cbed->scr =-1*numeric_limits<double>::max(); if (fhb.peek()!='\n') {hasScore=1;  fhb>>cbed->scr;}
    cbed->strand.clear();                     if (fhb.peek()!='\n') fhb>>cbed->strand;
    cbed->otro.clear();                       if (fhb.peek()!='\n') std::getline(fhb, cbed->otro);
    /* 
    printf("Load abed: %s %ld %ld %s",cbed->chr.c_str(),cbed->bgn,cbed->end, cbed->id.c_str());
    if (cbed->scr > numeric_limits<double>::min()) printf("  opt:%f %s %s", cbed->scr,cbed->strand.c_str(),cbed->otro.c_str());
    printf("\n");
   */ 
   if (ctree.find(cbed->chr) == ctree.end()) ctree[cbed->chr] = new IntervalTree();
    ctree[cbed->chr]->Insert(new myIntv(cbed->bgn,cbed->end,cbed->scr,cbed->id));
  }
  if (!fhb.eof()) cerr<<fn<<"  !!! not all parsed! Maybe croaked on score croaked\n";
  //fh.close();

  return ctree;
}

void print_bed(bed_s *b) {
  printf("%s\t%ld\t%ld\t%s",b->chr.c_str(),b->bgn-1,b->end, b->id.c_str());
  if (b->scr > -1*numeric_limits<double>::max()) printf("\t%f", b->scr); else if (dmyScore) printf("\t.");
  if (!b->strand.empty()) printf("\t%s", b->strand.c_str());
  if (!b->otro.empty()) printf("\t%s", b->otro.c_str());
}

int main ( int argc, char *argv[] ) {

  /* extensions
     1) can do multiple qrys
     2) % olap
     3) output each individual bedB or bedB scores
     4) output intersection
     5) 
     */

  int c;
  bool ocount = 0;
  bool oscore = 0;
  bool allscore = 0;
  bool allid = 0;
  while ( (c = getopt(argc, argv, "acsi")) != -1)   //c: count olaps, s: score overlaps, a: print scr of overlaps
    switch (c) {
      case 'c':  ocount = 1; break;
      case 's':  oscore = 1; break;
      case 'a':  allscore = 1; break;
      case 'i':  allid = 1; break;
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
      "-c  append the number of bedA lines olapping each bedB line [including 0]\n"<<
      "-s  append the sum of bedA line scores to each bedB line [including 0]\n"<<
      "-a  append all of the bedA line scores to each bedB line [comma separated inluding none]\n"<< 
      "-i  append all of the bedA line ids to each bedB line [comma separated inluding none]\n"<< 
      //"-b  append all of the bedA lines to each bedB line [comma separated]\n"<< 
      // just the intersection
      endl; 
    abort();
  } 
  string abed(argv[argc-2]); 
  string bbed(argv[argc-1]); 
  //printf("cnt:%d scr:%d abed:%s bbed:%s optind:%d %d\n",ocount,oscore,abed.c_str(),bbed.c_str(), optind,argc);
  if (abed == "-" && bbed=="-") {cout<<"both a and b bed can't be -\n";abort();}

  chromtree ctree = fn2IntervalTree(abed);

  ifstream fh(bbed.c_str()); if (!fh.is_open()) {cerr<<"Unable to open file: "<<bbed<<endl ;abort();}
 //fh.exceptions ( ifstream::failbit | ifstream::badbit );;

  bed_s *cbed = new bed_s();
  //  IntervalTree *tree;
  chromtree::iterator tree;



  bool print_all = (ocount||oscore||allscore||allid);
  while(fh>>cbed->chr>>cbed->bgn>>cbed->end) {

    cbed->bgn++;
    // optional fields   
    cbed->id.clear();                         if (fh.peek()!='\n') fh>>cbed->id;
    // score could be .
    cbed->scr =-1*numeric_limits<double>::max(); 
    if (fh.peek() !='\n'){hasScore=1; char d =fh.get();
      if (fh.peek() != '.'){ fh.putback(d); fh>>cbed->scr; }
      else                 { dmyScore=1;fh.get();}
    }
    cbed->strand.clear();                     if (fh.peek()!='\n') fh>>cbed->strand;
    cbed->otro.clear();                       if (fh.peek()!='\n') std::getline(fh, cbed->otro);



    if ( (tree = ctree.find(cbed->chr)) != ctree.end()) {
      // (Note: these aren't returned in order)
      TemplateStack<void *> *r = tree->second->Enumerate(cbed->bgn,cbed->end);

      // print bedB line
      if (print_all || r->Size() > 0) print_bed(cbed); 

      // append results
      if (ocount) { printf("\t%d",r->Size()); } 
      else if (oscore) {
        double score = 0;
        for (int j=0; j< r->Size(); j++) 
          if (((myIntv*)(*r)[j])->GetScr() > -1*numeric_limits<double>::max())
            score +=  ((myIntv*)(*r)[j])->GetScr();
        printf("\t%f",score);
      } else if (allscore) {
        int cnt = 0;
        for (int j=0; j< r->Size(); j++) 
          if (((myIntv*)(*r)[j])->GetScr() > -1*numeric_limits<double>::max()) {
            if (cnt++) printf(","); else printf("\t");
            printf("%f",((myIntv*)(*r)[j])->GetScr() );
          }
      }
     if (allid) {
        int cnt = 0;
        for (int j=0; j< r->Size(); j++)  {
          if ( !((myIntv*)(*r)[j])->GetId().empty() ) {
            if (cnt++) printf(","); else printf("\t");
            printf("%s",((myIntv*)(*r)[j])->GetId().c_str() );
          }
        }
      }
      if (print_all || r->Size() > 0) printf("\n");
     //r->Clear(delete);
     delete r;

    } else if (print_all) { print_bed(cbed); ocount ? printf("\t0\n") : printf("\n"); }
   //delete cbed->chr;delete cbed->id;delete cbed->strand;delete cbed->otro;
  }
  // can't close istream, not sure what to do
  //fh.close(); 
}
