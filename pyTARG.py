### pyTARG.py
#Author Sergio Bordel sergio_bordel@hotmail.com

#Permission is hereby granted, free of charge, to any person obtaining a copy 
#of this software and associated documentation files (the "Software"), to deal 
#in the Software without restriction, including without limitation the rights 
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell 
#copies of the Software, and to permit persons to whom the Software is furnished 
#to do so, subject to the following conditions:

#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

#import the cobra library
import cobra
from cobra import Model, Reaction, Metabolite



#Returns a flux distribution in mmol/h/g-biomass form a constrained genome scale metabolic model
def flux(model):
	sol=cobra.flux_analysis.optimize_minimal_flux(model)
	fluxes={}
	for r in sol.x_dict:
		fluxes[r]=0.027*sol.x_dict[r]
	return fluxes
	

def constrain(model,dic,lev,bound):
	mod=model
	for r in mod.reactions:
		if len(r.genes)>0:
			activo=1
			for gen in r.genes:
				if gen.id in dic:
					activo=0
					
			if activo==0:
				for ge in r.genes:
					
					if ge.id in dic and dic[ge.id]>lev:
						activo=1
			if activo==0:
				if r.upper_bound>0:
					r.upper_bound=bound
				if r.lower_bound<0:
					r.lower_bound=-bound
	return mod

#Constraints a model based on gene expression
def fullconstrain(model,dic,coef):
	mod=model
	for i in range(100):
		lev=10*range(100)[-i-1]
		bound=coef*lev
		
		mod=constrain(mod,dic,lev,bound)
	return mod

#Computes the effect of constraining one or several metabolic fluxes to 0.1 times their initial values
def block(model,targets):
	rea=[]
	sols=[]
	for r in model.reactions:
		if r.id in targets:
			rea.append(r)
		for g in r.genes:
			if g.id in targets:
				rea.append(r)
		for m in r.metabolites:
			if m.id in targets:
				rea.append(r)
	ubounds={}
	lbounds={}
	sol=cobra.flux_analysis.optimize_minimal_flux(model)
	reference=sol.f
	#print reference
	for r in rea:
		ubounds[r.id]=r.upper_bound
		lbounds[r.id]=r.lower_bound
	
	
	for r in rea:
		if sol.x_dict[r.id]>0:
			r.upper_bound=0.1*sol.x_dict[r.id]
		else:
			r.lower_bound=0.1*sol.x_dict[r.id]
	sol2=cobra.flux_analysis.optimize_minimal_flux(model)
	if reference>0.0000000000000000000001:
		result=sol2.f/reference
	else:
		result='not growth'
	
	
	return result
#Compares two models and returns a set of reactios that have a larger impact on the objective function of the first model while keeping as unaffected as possible the second
def personal(model,model2):
	sol1=cobra.flux_analysis.optimize_minimal_flux(model)
	sol2=cobra.flux_analysis.optimize_minimal_flux(model2)
	scores={}
	for r in model.reactions:
		if len(r.genes)>0:
			if abs(sol1.x_dict[r.id])>0.000001:
				score=abs(sol2.x_dict[r.id])-abs(sol1.x_dict[r.id])
				scores[r.id]=score
	sorted_vals=sorted(scores.items(), key=operator.itemgetter(1))
	ranked=[]
	for s in sorted_vals:
		ranked.append(s[0])
	refe=sol1.f
	refe1=sol1.f
	refe2=sol2.f
	resis=[]
	#find first reaction
	for r in ranked:
		if len(model.reactions.get_by_id(r).genes)>0:
			if abs(sol1.x_dict[r])>0.000001:
				up=model.reactions.get_by_id(r).upper_bound
				low=model.reactions.get_by_id(r).lower_bound
				if sol1.x_dict[r]>0:
					model.reactions.get_by_id(r).upper_bound=0.1*sol1.x_dict[r]
				else:
					model.reactions.get_by_id(r).lower_bound=0.1*sol1.x_dict[r]
				sol1=cobra.flux_analysis.optimize_minimal_flux(model)
				model.reactions.get_by_id(r).upper_bound=up
				model.reactions.get_by_id(r).lower_bound=low
				ratio=sol1.f/refe1
				
				if ratio<0.9:
					up=model2.reactions.get_by_id(r).upper_bound
					low=model2.reactions.get_by_id(r).lower_bound
					if sol2.x_dict[r]>0:
						model2.reactions.get_by_id(r).upper_bound=0.1*sol2.x_dict[r]
					else:
						model2.reactions.get_by_id(r).lower_bound=0.1*sol2.x_dict[r]
					sol2=cobra.flux_analysis.optimize_minimal_flux(model2)
					model2.reactions.get_by_id(r).upper_bound=up
					model2.reactions.get_by_id(r).lower_bound=low
					ratio2=sol2.f/refe2
					
					if ratio2>ratio+0.05:
						resis.append(r)
						
						refe1=sol1.f
						refe2=sol2.f
						print refe1
						print refe2
						if sol1.f/refe<0.5:
							break
						if sol1.x_dict[r]>0:
							model.reactions.get_by_id(r).upper_bound=0.1*sol1.x_dict[r]
						else:
							model.reactions.get_by_id(r).lower_bound=0.1*sol1.x_dict[r]
						if sol2.x_dict[r]>0:
							model2.reactions.get_by_id(r).upper_bound=0.1*sol2.x_dict[r]
						else:
							model2.reactions.get_by_id(r).lower_bound=0.1*sol2.x_dict[r]
				
	return resis


