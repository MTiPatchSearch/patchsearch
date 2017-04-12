#' @useDynLib patchsearch
#' @importFrom Rcpp evalCpp

patchsearch_MINCLIQUE=6
patchsearch_MAXCLIQUE=10
patchsearch_BCMIN=0.
patchsearch_MINDIST=0.2
patchsearch_MAXDIST=3
patchsearch_SEQHOM=0
patchsearch_MERGECLIQUE=1
patchsearch_NBNEI=4
patchsearch_LOCDEGREE=4
patchsearch_LOCSIZE=6
patchsearch_ATYPES=c("O","N","S","CA","Ca","C")

#' extract patch from a pdb file
#'
#' @param chains
#' @param ligchain
#' @param ligresno
#' @param cutoff
#' @export

patch_extract=function(pdb, ligpdb, prot.inds, lig.inds, cutoff) {
    cat("patch_extract...\n")
    if (length(lig.inds$atom)==0) {
       message ("find no ligand..")
       return(NULL)
    }
    bs = bio3d::binding.site(pdb,ligpdb,
		      a.inds=prot.inds, b.inds=lig.inds,
		      cutoff=cutoff)
    #ind=pdb$atom$eleno[bs$ind$atom]
    #if (add.cycles) ind=arom.add(pdb, ind)
    #cat("patch_extract::add.calpha",add.calpha,"\n")
    #if (add.calpha) ind=calpha.add(pdb, ind)
    patch.pdb=bio3d::trim.pdb(pdb, resno=bs$resno)
    return(patch.pdb)
}


chain_surf_type=function(surf.pdb, pdb, chain) {
    surf.pdb=bio3d::trim.pdb(surf.pdb, chain=chain)
    if (any(is.na(surf.pdb$atom$elesy))) {
       surf.pdb$atom$elesy=substring(surf.pdb$atom$elety,1,1) # si pas de elesy
    }
    surf.pdb$atom$elesy[surf.pdb$calpha]="CA"
    resno=unique(surf.pdb$atom$resno)
    ca=bio3d::trim.pdb(surf.pdb, resno=resno, chain=chain, string="calpha")
    for (r in 1:length(ca$atom$eleno))  {
      resid=ca$atom$resid[r]
      ch=ca$atom$chain[r]
      resno=ca$atom$resno[r]
      if (resid=="PHE") {
        ind=bio3d::atom.select(pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
       #pdb$atom[ind$atom, ]
       	if (length(ind$xyz)/3>=6) {
	      center=bio3d::com(pdb, inds=ind, use.mass=FALSE)
              inds=bio3d::atom.select(surf.pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
	      if (length(inds$atom)>0) {
     	          ind1=inds$atom[1]
       	      	  surf.pdb$atom$x[ind1]=center[1]
       	      	  surf.pdb$atom$y[ind1]=center[2]
       	      	  surf.pdb$atom$z[ind1]=center[3]
       	      	  surf.pdb$atom$elesy[ind1]="Ca"
     	      	  if (length(inds$atom)>1) surf.pdb$atom=surf.pdb$atom[-inds$atom[-1], ]
	      }
       	 }
      }
      if (resid=="TYR") {
       	   ind=bio3d::atom.select(pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
      	   if (length(ind$xyz)/3>=6) {
       	      center=bio3d::com(pdb, inds=ind, use.mass=FALSE)
              inds=bio3d::atom.select(surf.pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE1","CE2","CZ"))
      	      ind1=inds$atom[1]
	      if (length(inds$atom)>0) {
       	           surf.pdb$atom$x[ind1]=center[1]
      	      	   surf.pdb$atom$y[ind1]=center[2]
       	      	   surf.pdb$atom$z[ind1]=center[3]
       	      	   surf.pdb$atom$elesy[ind1]="Ca"
       	      	   if (length(inds$atom)>1) surf.pdb$atom=surf.pdb$atom[-inds$atom[-1],]
 	      }
      	   }
      }
      if (resid=="TRP") {
       	   # replacement of first cycle
       	   ind=bio3d::atom.select(pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE2","NE1"))
      	   if (length(ind$xyz)/3>=5) {
              center=bio3d::com(pdb, inds=ind, use.mass=FALSE)
       	      inds=bio3d::atom.select(surf.pdb, resno=resno, chain=ch, elety=c("CG","CD1","CD2","CE2","NE1"))
 	      if (length(inds$atom)>0) {
     	      	 ind1=inds$atom[1]
		 surf.pdb$atom$x[ind1]=center[1]
       	      	 surf.pdb$atom$y[ind1]=center[2]
       	      	 surf.pdb$atom$z[ind1]=center[3]
       	      	 surf.pdb$atom$elesy[ind1]="Ca"
       	      	 if (length(inds$atom)>1) surf.pdb$atom=surf.pdb$atom[-inds$atom[-1], ]
	      }
	   }
       	   # replacement of second cycle
       	   ind=bio3d::atom.select(pdb, resno=resno, chain=ch, elety=c("CD2","CE2","CE3","CZ2","CH2","CZ3"))
      	   if (length(ind$xyz)/3>=6) {
       	      center=bio3d::com(pdb, inds=ind, use.mass=FALSE)
       	      inds=bio3d::atom.select(surf.pdb, resno=resno, chain=ch, elety=c("CD2","CE2","CE3","CZ2","CH2","CZ3"))
  	      if (length(inds$atom)>0) {
     	      	 ind1=inds$atom[1]
       	      	 surf.pdb$atom$x[ind1]=center[1]
       	      	 surf.pdb$atom$y[ind1]=center[2]
       	      	 surf.pdb$atom$z[ind1]=center[3]
       	      	 surf.pdb$atom$elesy[ind1]="Ca"
      	      	 if (length(inds$atom)>1) surf.pdb$atom=surf.pdb$atom[-inds$atom[-1], ]
	      }
     	   }
       }
   }
   surf.pdb$xyz=bio3d::as.xyz(as.vector(t(surf.pdb$atom[,c("x","y","z")])))
   surf.pdb
}

#' type atoms on a surface in a pdb file
#'
#' @param pdb
#' @export
surf_type=function(surf.pdb, pdb, chains) {
	out.pdb=chain_surf_type(surf.pdb, pdb, chains[1])
	if (length(chains)>=2) 
	    for (k in 2:length(chains)) {
	    	tmp.pdb=chain_surf_type(surf.pdb, pdb, chains[k])
		    out.pdb=suppressWarnings(bio3d::cat.pdb(out.pdb, tmp.pdb, rechain=FALSE, renumber=FALSE))
		}
	out.pdb
}


naccess=function(pdb) {
    id=tempfile(tmpdir=getwd())
    pdbfile=paste(id, ".pdb", sep="")
    bio3d::write.pdb(pdb, file=pdbfile)
    command=paste("naccess -a", pdbfile)
    system(command)
    A=suppressWarnings(bio3d::read.pdb(paste(id,".asa",sep="")))
    acc=A$atom[,"o"]
    names(acc)=paste(A$atom$elety,A$atom$resid, sep="_")
    acc=acc/patchsearch_atom_ACCS[names(acc)]
    acc[is.na(acc)]=0
    unlink("*.asa")
    unlink("*.log")
    unlink(pdbfile)
    acc
}

# ind: eleno=number of selected atoms (surface or patch atoms)
calpha.add=function(pdb, ind) {
    sel=bio3d::atom.select(pdb, eleno=ind)
    selres=unique(pdb$atom$resno[sel$atom])
    sel=bio3d::atom.select(pdb, string="calpha", resno=selres)
    ind=union(ind,pdb$atom$eleno[sel$atom])
    o=order(ind)
    ind[o]
}

naccess_chain_surf_extract=function(pdb, chain, minacc, add.calpha) {
    cat("naccess chain_surf_extract=", chain, "\n")
    ind1=bio3d::atom.select(pdb,  string="protein", chain=chain)
    ind2=bio3d::atom.select(pdb,  string="noh", chain=chain)
    ind=bio3d::combine.select(ind1, ind2)
    new.pdb=bio3d::trim.pdb(pdb, inds=ind)
    acc=naccess(new.pdb)
    ind=new.pdb$atom$eleno[acc>=minacc] # seleted indices of atom (eleno)
    if (add.calpha) ind=calpha.add(new.pdb, ind)
    selatom=ind
    surf.pdb=bio3d::trim.pdb(new.pdb,eleno=selatom)
    return(surf.pdb)
}

dssp_chain_surf_extract=function(pdb, chain, minacc) {
    cat("dssp chain_surf_extract=", chain, "\n")
    new.pdb=bio3d::trim.pdb(pdb, string="protein", chain=chain)
    result=bio3d::dssp(new.pdb)
    acc=result$acc
    ca.pdb=bio3d::trim.pdb(pdb, string="calpha", chain=chain)
    names(acc)=tolower(ca.pdb$atom$resid)
    acc=acc/patchsearch_ACCS[names(acc)]
    selres=ca.pdb$atom$resno[acc>=minacc]
    surf.pdb=bio3d::trim.pdb(new.pdb,resno=selres)
    return(surf.pdb)
}

chain_surf_extract=function(pdb, chain, minacc, add.calpha, method) {
    surf.pdb=NULL
    if (method=="all") surf.pdb=pdb
    if (method=="dssp") surf.pdb=dssp_chain_surf_extract(pdb, chain, minacc)
    if (method=="naccess") surf.pdb=naccess_chain_surf_extract(pdb, chain, minacc, add.calpha)
    surf.pdb
}

#' extract a surface from a pdb file given a list of chains
#'
#' @param pdb
#' @param chains
#' @param minacc
#' @export
surf_extract=function(pdb, chains, minacc=MINACC,  add.calpha=TRUE, method="naccess") {
    cat("surf_extract...\n")
    surf.pdb=chain_surf_extract(pdb, chains[1], minacc, add.calpha, method)
    if (length(chains)>1) {
	for (i in 2:length(chains)) {
	    tmp.pdb=chain_surf_extract(pdb, chains[i], minacc, add.calpha, method)
	    surf.pdb=bio3d::cat.pdb(surf.pdb, tmp.pdb, renumber=FALSE, rechain=FALSE)
	    # renumber=FALSE doesn't work when chain break
	 }
    }
    return(surf.pdb)
}


#' Compute BC score between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
BCscore<-function(X,Y) {
    X=scale(X,scale=FALSE)
    Y=scale(Y,scale=FALSE)
    if (det(t(X)%*%X)<1 | det(t(Y)%*%Y)<1) return(0)
    det(t(X)%*%Y)/(sqrt(det(t(X)%*%X))*sqrt(det(t(Y)%*%Y)))
}

#' Compute rmsd between X and Y
#'
#' @param X matrix Nx3
#' @param Y matrix Nx3
#' @export
patchsearch_rmsd=function(X,Y) {
   r=bio3d::rmsd(as.vector(t(X)),as.vector(t(Y)),fit=TRUE)
   dev=devsub(X,as.integer(0:(nrow(X)-1)),Y,as.integer(0:(nrow(Y)-1)))
   #cat("dev=", sqrt(mean(dev**2)),"\n")
   #return (list(rmsd=r,dev=sqrt(abs(r**2-mean(dev)**2))))
   return (list(rmsd=r,dev=0))
}

#' read PDB file containing an extracted protein surface or patch
#'
#' @param file  pdb file
#' @export 
read_surface <- function(pdbfile) {
	pdb=suppressWarnings(bio3d::read.pdb(pdbfile))
	pdb=bio3d::trim.pdb(pdb, string="protein", resid=patchsearch_AA)
	x=pdb$atom$x
	y=pdb$atom$y
	z=pdb$atom$z
	type=gsub(" ","",pdb$atom$elesy)
	#n=length(type)+5*sum(type=="Ca")
	n=length(type)
	resname=patchsearch_codes[pdb$atom$resid]
	ch=pdb$atom$chain
	ch[is.na(ch)]=" "
    list(coord=cbind(x,y,z),type=type,num=pdb$atom$eleno,resnum=pdb$atom$resno,resname=resname, natoms=n, aname=pdb$atom$elety, chains=ch)
}


extract_cliques<-function(gp, N, M) {
        #Cl=largest.cliques(gp)
	size=igraph::clique.number(gp)
	cat("extract_cliques::size=",size,"\n")
        Cl=igraph::maximal.cliques(gp,min=size-1)
	#Cl=cliques(gp, min=6,max=6)
	if (length(Cl)==0) return(list(clique=NULL,size=0))
	##
	for (i in 1:length(Cl)) Cl[[i]]=as.integer(names(Cl[[i]]))
	##
	list(clique=Cl, size=size-1)
}	

union_clique <- function(C1,C2,N) {
  I1=(C1-1)%%N+1
  I2=(C2-1)%%N+1
  J1=(C1-1)%/%N+1
  J2=(C2-1)%/%N+1
  I=which(!I2%in%I1 & !J2%in%J1)
  base::union(C1,C2[I])
}

merge_cliques<-function(X,Y,clusters,N,common=patchsearch_MERGECLIQUE, bcmin=patchsearch_BCMIN) {
   stop=FALSE
   while (!stop) {
         stop=TRUE
   	 n=length(clusters)
	 i=1
  	 while (i<= n-1) {
	   j=i+1
	   while (j<= n) {
	     if (length(base::intersect(clusters[[i]],clusters[[j]]))>=common) {
	     	#clusters[[i]]=union(clusters[[i]],clusters[[j]]) # may create forks
	     	C=union_clique(clusters[[i]],clusters[[j]],N)
   		I=(C-1)%%N+1 # target indices
   		J=(C-1)%/%N+1# query indices
   		bc=BCscore(X[I,], Y[J,])
		if (bc>=bcmin) {
		   clusters[[i]]=C
		   clusters[[j]]=c()
		   n=n-1
		   stop=FALSE
		}
	     }
	     j=j+1
     	    }
	    i=i+1
	  }
   }
   return(clusters)
}

merge_cliques2<-function(X,Y,clusters,N,common=patchsearch_MERGECLIQUE) {
   stop=FALSE
   while (!stop) {
         stop=TRUE
   	 n=length(clusters)
	 i=1
  	 while (i<= n-1) {
	   j=i+1
	   while (j<= n) {
	     if (length(base::intersect(clusters[[i]],clusters[[j]]))>=common) {
	     	#clusters[[i]]=union(clusters[[i]],clusters[[j]]) # may create forks
	     	C=union_clique(clusters[[i]],clusters[[j]],N)
 		clusters[[i]]=C
		clusters[[j]]=c()
		n=n-1
		stop=FALSE
	     }
	     j=j+1
     	    }
	    i=i+1
	  }
   }
   return(clusters)
}

build_graph=function(types, mindist, X, XProp, Xseq, Y, YProp, Yseq, S, SEQHOM, locdegree=patchsearch_LOCDEGREE, locsize=patchsearch_LOCSIZE) {
	N=nrow(X)
	motif=paste(types,collapse="|")
	J=grep(motif, YProp)
	I=grep(motif, XProp)
	Mmotif=length(J)
	Nmotif=length(I)
	V=vertexc(XProp[I], Xseq[I], YProp[J], Yseq[J], S, SEQHOM)
	nV=nrow(V)
	locdegrees=integer(nV)
	if (all(V==0)) return(NULL)
	V=cbind(J[V[,1]],I[V[,2]]) # atom ids 1..M, 1..N 
	E=buildGraphc(X, Y, V, locdegrees, as.double(mindist), locsize)
	if (is.null(E)) return(NULL)
	P=cbind((E[,3]-1)*N+E[,1],(E[,4]-1)*N+E[,2])
	mode(P)="character" #vertex label and not id to avoid huge graphs with non connected nodes
        gp=igraph::graph.edgelist(P, directed=FALSE)
	#cat("graph:vertices:", igraph::vcount(gp), "graph.edges:", igraph::ecount(gp), length(locdegree),"N=",N,"\n")
	if (locdegree>0) {
	   delV=which(locdegrees<=locdegree)
	   delE=as.character((V[delV,1]-1)*N+V[delV,2])
	   gp=igraph::delete_vertices(gp, delE)
	   #cat("graph:vertices:", igraph::vcount(gp), "graph.edges:", igraph::ecount(gp),"\n")
	}
	return(gp)
}

clique_patchsearch=function(g, X, Y, minclique=patchsearch_MINCLIQUE, bcmin=patchsearch_BCMIN) {
	N=nrow(X);M=nrow(Y)
	cl=extract_cliques(g,N,M)
	##
	#cat("clique_patchsearch::cl=",cl$size,"\n")
	if (cl$size<minclique) return(NULL)
	nbclique=length(cl$clique)
	size=bc=double(nbclique)
	for (i in 1:nbclique) {
	    C=cl$clique[[i]]
	    size[i]=length(C)
	    I=(C-1)%%N+1 # target indices
	    J=(C-1)%/%N+1# query indices
	    bc[i]=BCscore(X[I,], Y[J,])
	}
	#cat("clique_patchsearch::bc=",sum(bc>=bcmin),"\n")
	cl$clique=cl$clique[bc>=bcmin]
	#
	if (length(cl$clique)==0) return(NULL)
	bc=bc[bc>=bcmin]
	o=order(bc,decreasing=TRUE)
	cl$bc=bc[o]
	cl$clique=cl$clique[o]
	cl$size=size[o]
	return(cl)
}



#' search for matchings between query pdb and a surface pdb files
#'
#' @param query
#' @param target
#' @param 
#' @export
psearch=function(query, target, mindist=patchsearch_MINDIST, maxdist=patchsearch_MAXDIST, minclique=patchsearch_MINCLIQUE, mergeclique=patchsearch_MERGECLIQUE, bcmin= patchsearch_BCMIN, seqhom=patchsearch_SEQHOM, enrich=TRUE, merge=TRUE, ptypes=patchsearch_ATYPES,locdegree=patchsearch_LOCDEGREE, locsize=patchsearch_LOCSIZE) {
    #
    # cat("readPDBs\n")
    #
    X=target$coord
    Y=query$coord
    XProp=target$type
    YProp=query$type
    Xseq=target$resname
    Yseq=query$resname
    NN=target$natoms
    MM=query$natoms
    N=nrow(X)
    M=nrow(Y)
    #
    # construction du graphe produit
    #
    #
    # recherche des cliques
    #
    cliques=NULL
    k=length(ptypes)
    atypes=c(ptypes, setdiff(patchsearch_ATYPES,ptypes)) # first ptypes then others
    while (k<=length(atypes) & is.null(cliques)) {
	  cat("graph=",atypes[1:k],"\n")
    	  graph=build_graph(atypes[1:k], mindist, X, XProp, Xseq, Y, YProp, Yseq, patchsearch_S, seqhom, locdegree, locsize)
	  cat("clique search...\n")
    	  cliques=clique_patchsearch(graph, X, Y, minclique, bcmin)
	  k=k+1
    }
    k=k-1
    motif=paste(atypes[1:k],collapse="|")
    J=grep(motif, YProp)
    atom_nb=length(J)
    cat("number of cliques=",length(cliques$clique),"size max=", max(cliques$size), "total number of atoms=",atom_nb,"\n")
    imax=which.max(cliques$bc*cliques$size/atom_nb)
    rbc0=cliques$bc[imax]*cliques$size[imax]/atom_nb
    alen0=cliques$size[imax]
    cat("best clique:","size=", alen0, "rbc=", rbc0, "atom_nb=", atom_nb,"\n")
    if (is.null(cliques)) {
	message("no valid clique")
	return(NULL)
    }
    cliques=cliques$clique
    #
    # merging cliques
    #
    #clusters=merge_cliques(X,Y,cliques,N, mergeclique,-1)
    clusters=merge_cliques2(X,Y,cliques,N, mergeclique)
    #clusters=cliques
    nbclique=length(clusters)
    cat("nb of clique clusters=",nbclique,"\n")

    bc=rbc=double(nbclique)
    for (i in 1:nbclique) {
    	C=clusters[[i]]
   	size=length(C)
   	I=(C-1)%%N+1 # target indices
   	J=(C-1)%/%N+1# query indices
   	bc[i]=BCscore(X[I,], Y[J,])
	rbc[i]=size/atom_nb*bc[i]
    }
    clusters=clusters[bc>=bcmin]
    rbc=rbc[bc>=bcmin]
    bc=bc[bc>=bcmin]
    o=order(rbc,decreasing=TRUE)
    clusters=clusters[o]
    bc=bc[o]
    rbc=rbc[o]
    nbclique=length(clusters)
    if (nbclique==0) {
	message("no valid clique after merging")
	return(NULL)
    }
    cat("nb of selected clique clusters=",nbclique,"\n")
#
# enriching the MAXCLIQUES best valid cliques 
#
if (enrich==TRUE) {
   V=vertexc(XProp, Xseq, YProp, Yseq, patchsearch_S, -10)
   for (ic in 1:nbclique) {
       C=clusters[[ic]]
       nbefore= length(C)
    #
    # calcul de K=ensemble de links à tester pour enrichir clique choisie
    # 1. ensemble des voisins de C
    #
       cat("enrichment...", length(C),"\n")
       K=getNeighborsc(as.integer(C),V,X,Y,as.double(maxdist),as.integer(patchsearch_NBNEI))
       cat("enrichment first pass...\n")
       Cplus=enrichloopc(as.integer(C), as.integer(K), X, Y, as.double(maxdist))
       cat("enrichment second pass...", length(Cplus),"\n")
       C=Cplus
       Cplus=enrichloopc(as.integer(C), as.integer(K), X, Y, 100.)
       cat("enrichment end...", length(Cplus),"\n")
       clusters[[ic]]=Cplus
    }
}
nclique=length(clusters)
cat("nb of enriched quasi-cliques=",length(clusters),"\n")
if (merge==TRUE) {
stop=FALSE
    while(!stop) {
    	clusters=merge_cliques2(X,Y,clusters,N, 6)
    	cat("nb clustered enriched quasi-cliques=",nbclique,"\n")
    	if (nbclique==length(clusters)) stop=TRUE
    	nbclique=length(clusters)
    }
}
#
# returns matched atoms as indices of query (Iquery) and target (Itarget)
#
   Itarget=list()
   Iquery=list()
   rbc=double(nbclique)
   rms=double(nbclique)
   maxr=double(nbclique)
   alen=integer(nbclique)
   for (ic in 1:nbclique) {
       I=(clusters[[ic]]-1)%%N+1 # target indices
       J=(clusters[[ic]]-1)%/%N+1# query indices
       Itarget[[ic]]=I
       Iquery[[ic]]=J
       alen[ic]=length(J)
       rbc[ic]=alen[ic]/MM*BCscore(X[I,], Y[J,])
       r=patchsearch_rmsd(X[I,],Y[J,])
       rms[ic]=r$rmsd
       maxr[ic]=r$dev
   }
   return(list(target=Itarget, query=Iquery, alen=alen, rmsd=rms, rbc=rbc, maxr=maxr, alen0=alen0, rbc0=rbc0))
}
