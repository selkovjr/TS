/*
 * Copyright (C) 2010 Life Technologies Corporation. All rights reserved.
 */

/*
 * prune.hpp
 *
 *  Created on: Mar 23, 2010
 *      Author: kennedcj
 */

#ifndef PRUNE_HPP_
#define PRUNE_HPP_

#include "node.hpp"
#include <fstream>

namespace relation {
  template<class Label, class Weight=int>
  class None {
 	  public:
			enum Kind {NONE, PRINT, ERASE};

 	    None(void) : p_kind(NONE) {std::cout << "None()" << std::endl;};

 		  virtual bool operator()(Node<Label, Weight> const& node) {return false;}
 		  virtual bool operator()(typename Node<Label, Weight>::iterator edge) {return false;}

 		  virtual std::string getKind(void) const {return "None::";};

 		  Kind kind(void) const {return p_kind;}

 		  virtual ~None(void) {};

 	  protected:
 		  Kind p_kind;
  };

  template<class Label, class Weight=int>
  class Print : public None<Label, Weight> {
    public:
  	  // Print(void) {/*fp = &std::cout;*/};
  	  enum Format {SIF, DOT};

  	  Print(Format f=SIF, const char* file=NULL) : format(f) {
  	  	this->p_kind = None<Label, Weight>::PRINT;

  	  	if(file) {
  	  		fout = new std::ofstream(file);
  	  		assert(fout->is_open());
  	  		fp = fout;
  	  	} else {
  	  		fout = NULL;
  	      fp = &std::cout;
  	      std::cout << "fp = std::cout ~" << std::endl;
  	  	}

  	  	switch(format) {
  	  	  case DOT: {
  	  	  	*fp << "digraph G {" << std::endl;
  	  	  } break;
  	  	}

  	  }



  	  // using None<Label, Weight>::operator();
  	  virtual bool operator()(Node<Label, Weight> const& node) {
  	  	typename Node<Label, Weight>::iterator edge = node.begin(), last = node.end();

  	  	switch(format) {
  	  	  case SIF: {
  	  	    while(edge != last) {
  	  	  	  *fp << node.label() << " " << (*edge).weight() << " " << (*edge).target().label() << std::endl;
  	  	  	  ++edge;
  	  	  	}
  	  	  } break;
  	  	  case DOT: {
  	  	  	while(edge != last) {
  	  	  		*fp << node.label() << " -> " << (*edge).target().label() << " [label=" << (*edge).weight() << "];" << std::endl;
  	  	  		++edge;
  	  	  	}
  	  	  } break;
  	  	}

  	  	return false;
  	  }

  	  using None<Label, Weight>::operator();
  	  //bool operator()(typename Node<Label, Weight>::iterator edge) {return false;}


  	  // using None<Label, Weight>::getKind;
  	  std::string getKind(void) const {return "Print::";};


  	  //bool operator()(typename Node<Label, Weight>::iterator edge) {
  	  //	return false;
  	  //}

  	  ~Print(void) {
  	  	if(fout) {
    	  	switch(format) {
    	  	  case DOT: {
    	  	  	*fp << "}" << std::endl;
    	  	  } break;
    	  	}

  	  		fout->close();
  	  		// delete fout; // TODO: is this a memory leak? Deleting causes seg fault.
  	  	}
  	  }

  	private:
  	  Format format;
  	  std::ostream *fp;
  	  std::ofstream *fout;
  };
}

#endif /* PRUNE_HPP_ */
