#ifndef DUNE_NONLINEARPARABOLIC_HH
#define DUNE_NONLINEARPARABOLIC_HH

#include<dune/disc/shapefunctions/lagrangeshapefunctions.hh>
#include<dune/disc/functions/p1function.hh>
#include<dune/disc/operators/p1operator.hh>
#include"dumux/nonlinear/nonlinearmodel.hh"

namespace Dune
{
  template<class G, class RT, class ProblemType, class LocalJacobian, 
            class FunctionType, class OperatorAssembler>
  class NonlinearParabolic 
  : public NonlinearModel<G, RT, ProblemType, LocalJacobian, FunctionType, OperatorAssembler> 
  {
  public:	
	typedef NonlinearModel<G, RT, ProblemType, LocalJacobian, 
	                          FunctionType, OperatorAssembler> NonlinearModel;
	
	NonlinearParabolic(const G& g, ProblemType& prob)
	: NonlinearModel(g, prob), uOldTimeStep(g)
	{ }
	
	NonlinearParabolic(const G& g, ProblemType& prob, int level)
	: NonlinearModel(g, prob, level), uOldTimeStep(g, level)
	{ 	}
	
	virtual void initial() = 0;
	
	virtual void update(double& dt) = 0;
	
	virtual void solve() = 0;
	
	FunctionType uOldTimeStep;
  };


  
  
  
  template<class G, class RT, class ProblemType, class LocalJac, int m=1>
  class LeafP1NonlinearParabolic : public NonlinearParabolic<G, RT, ProblemType, LocalJac, 
                                        LeafP1Function<G, RT, m>, LeafP1OperatorAssembler<G, RT, m> >
  {
  public:
	  // define the function type:
	  typedef LeafP1Function<G, RT> FunctionType;

	  // define the operator assembler type:
	  typedef LeafP1OperatorAssembler<G, RT, m> OperatorAssembler;

	  typedef NonlinearParabolic<G, RT, ProblemType, LocalJac, 
	                          FunctionType, OperatorAssembler> NonlinearParabolic;
	  
	  typedef LeafP1NonlinearParabolic<G, RT, ProblemType, LocalJac, m> ThisType;

	  typedef LocalJac LocalJacobian;
	  
	  // mapper: one data element per vertex
	  template<int dim>
	  struct P1Layout
	  {
		  bool contains (Dune::GeometryType gt)
		  {
			  return gt.dim() == 0;
		  }
	  }; 

	  typedef typename G::Traits::LeafIndexSet IS;
	  typedef MultipleCodimMultipleGeomTypeMapper<G,IS,P1Layout> VertexMapper;

	  LeafP1NonlinearParabolic (const G& g, ProblemType& prob) 
	  : NonlinearParabolic(g, prob), grid(g), vertexmapper(g, g.leafIndexSet())
	  { }
	  
	  virtual void initial() 
	  {
		  typedef typename G::Traits::template Codim<0>::Entity Entity;
		  typedef typename G::ctype DT;
		  typedef typename IS::template Codim<0>::template Partition<All_Partition>::Iterator Iterator;
		  enum{dim = G::dimension};
		  enum{dimworld = G::dimensionworld};
		  
		  const IS& indexset(grid.leafIndexSet());
		  
		  // iterate through leaf grid an evaluate c0 at cell center
		  Iterator eendit = indexset.template end<0, All_Partition>();
		  for (Iterator it = indexset.template begin<0, All_Partition>(); it != eendit; ++it)
		  {
			  // get geometry type
			  Dune::GeometryType gt = it->geometry().type();

		      const typename Dune::LagrangeShapeFunctionSetContainer<DT,RT,dim>::value_type& 
		      	sfs=Dune::LagrangeShapeFunctions<DT,RT,dim>::general(gt, 1);
		      int size = sfs.size();

		      for (int i = 0; i < size; i++) {
		    	  // get cell center in reference element
		    	  const Dune::FieldVector<DT,dim>& 
					  local = sfs[i].position();

		    	  // get global coordinate of cell center
		    	  Dune::FieldVector<DT,dimworld> global = it->geometry().global(local);

		    	  int globalId = vertexmapper.template map<dim>(*it, sfs[i].entity());
		    	  
		    	  // initialize cell concentration
		    	  (*(this->u))[globalId] = this->problem.initial(global, *it, local);
		      }
		  }

		  *(this->uOldTimeStep) = *(this->u);
		  return;
	  }

		void vtkout (const char* name, int k) const 
		{
			VTKWriter<G, typename G::template Codim<0>::LeafIndexSet> 
				vtkwriter(this->grid, this->grid.leafIndexSet());
			char fname[128];	
			sprintf(fname,"%s-%05d",name,k);
			vtkwriter.addVertexData(*(this->u),"total pressure p~");
			vtkwriter.write(fname, VTKOptions::ascii);		
		}

  protected:
	  const G& grid;
	  VertexMapper vertexmapper;
  };

}
#endif
