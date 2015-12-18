#include "exportClasses.h"
#include "TObject.h"
#include "reader.h"
//#include "compact-7.3.x"

//##########################
//###   RawObjects
//##########################
ClassImp(NVtxTrack);
ClassImp(NSCVertex);
ClassImp(NTrak);
ClassImp(NCluster);
//##########################
//###   Corrected Objects
//##########################
ClassImp(NPhysicsCluster);
ClassImp(NPhysicsTrack);
//##########################
//###   Physics Objects
//##########################
ClassImp(NRecoParticle);
ClassImp(NMCParticle);
//##########################
//###   Database
//##########################
ClassImp(NSuperTimeOffset);
ClassImp(NAbcog_params);
ClassImp(NGeom);
ClassImp(NDCH);
ClassImp(Nxyz);
ClassImp(NDETStatus);

//##########################
//###   Top nodes
//##########################
ClassImp(ROOTRawEvent);
ClassImp(ROOTCorrectedEvent);
ClassImp(ROOTBurst);
ClassImp(ROOTFileHeader);
ClassImp(ROOTPhysicsEvent);
ClassImp(ROOTMCEvent);

NVtxTrack::NVtxTrack(vtxtracks &ref):
		iTrack(ref.iTrack), bdxdz(ref.bdxdz), bdydz(ref.bdydz)
{
};

NSuperTimeOffset& NSuperTimeOffset::operator=(superTimeOffset &ref){
	Version = ref.Version;
	Tag = ref.Tag;
	Aks = ref.Aks;
	Kab = ref.Kab;
	Nmv = ref.Nmv;
	Akl = ref.Akl;
	Hod = ref.Hod;
	Nho = ref.Nho;
	Dch = ref.Dch;
	Lkr = ref.Lkr;
	Hac = ref.Hac;
	Muv = ref.Muv;
	LkrTag = ref.LkrTag;
	LkrNhod = ref.LkrNhod;
	LkrAkl = ref.LkrAkl;
	LkrHac = ref.LkrHac;
	KabPlus = ref.KabPlus;
	KabMinus = ref.KabMinus;

	return *this;
};

NAbcog_params& NAbcog_params::operator=(void *r)
{
	abcog_params_t *ref = (abcog_params_t*)r;

	alpha = ref->alpha;
	alpha_coeff = ref->alpha_coeff;
	beta = ref->beta;
	beta_coeff = ref->beta_coeff;
	mkp = ref->mkp;
	mkperr = ref->mkperr;
	mkn = ref->mkn;
	mknerr = ref->mknerr;
	cogX1p = ref->cogX1p;
	cogY1p = ref->cogY1p;
	cogX1n = ref->cogX1n;
	cogY1n = ref->cogY1n;
	cogX4p = ref->cogX4p;
	cogY4p = ref->cogY4p;
	cogX4n = ref->cogX4n;
	cogY4n = ref->cogY4n;
	status = ref->status;
	pkp = ref->pkp;
	pkdxdzp = ref->pkdxdzp;
	pkdydzp = ref->pkdydzp;
	pkxoffp = ref->pkxoffp;
	pkyoffp = ref->pkyoffp;
	pkm = ref->pkm;
	pkdxdzm = ref->pkdxdzm;
	pkdydzm = ref->pkdydzm;
	pkxoffm = ref->pkxoffm;
	pkyoffm = ref->pkyoffm;

	return *this;
};

ROOTRawEvent& ROOTRawEvent::operator=(superCmpEvent *ref){
	clear();
	Nvtx = ref->Nvtx;
	Ncluster = ref->Ncluster;
	timeStamp = ref->timeStamp;
	trigWord = ref->trigWord;

	for(unsigned int i=0; i<Nvtx; ++i){
		vtx.push_back(NSCVertex(ref->vtx[i]));
	}

	DETStatus = &ref->DETstatus[0];
	return *this;
}

ROOTBurst& ROOTBurst::operator=(superBurst *ref){
	nrun = ref->nrun;
	time = ref->time;
	tOffst = ref->tOffst;
	return *this;
}


NSCVertex::NSCVertex(SCvertex &ref):
	Nvtxtrack(ref.Nvtxtrack),
	charge(ref.charge),
	cda(ref.cda),
	chi2(ref.chi2),
	time(0),
	position(ref.x, ref.y, ref.z)
{
	for(unsigned int i=0; i<Nvtxtrack; ++i){
		vtxtrack.push_back(NVtxTrack(ref.vtxtrack[i]));
	}

}

NTrak::NTrak(trak &ref):
	p(ref.p),
	q(ref.q),
	quality(ref.quality),
	chi2(ref.chi2),
	by(ref.by),
	bx(ref.bx),
	bdxdz(ref.bdxdz),
	bdydz(ref.bdydz),
	vtxID(-1),
	time(ref.time),
	dDeadCell(ref.dDeadCell),
	bDetPos(ref.bx, ref.by, 0),
	aDetPos(ref.x, ref.y, 0),
	bMomentum(TVector3(ref.bdxdz, ref.bdydz, 1).Unit()),
	aMomentum(TVector3(ref.dxdz, ref.dydz, 1).Unit())
{

}

NCluster::NCluster(cluster &ref):
	time(ref.time),
	dDeadCell(ref.dDeadCell),
	E(ref.energy),
	lkr_acc(false),
	position(ref.x, ref.y, 0)
{
}


NGeom& NGeom::operator=(GeomCompact *ref){
	Dch[0].PosChamber.SetXYZ(ref->Dch[0].PosChamber.x, ref->Dch[0].PosChamber.y, ref->Dch[0].PosChamber.z);
	Dch[1].PosChamber.SetXYZ(ref->Dch[1].PosChamber.x, ref->Dch[1].PosChamber.y, ref->Dch[1].PosChamber.z);
	Dch[2].PosChamber.SetXYZ(ref->Dch[2].PosChamber.x, ref->Dch[2].PosChamber.y, ref->Dch[2].PosChamber.z);
	Dch[3].PosChamber.SetXYZ(ref->Dch[3].PosChamber.x, ref->Dch[3].PosChamber.y, ref->Dch[3].PosChamber.z);

	Lkr.SetXYZ(ref->Lkr.x, ref->Lkr.y, ref->Lkr.z);

	return *this;
}

NDETStatus& NDETStatus::operator =(DETstatus* ref)
{
	//TAG = ref->TAG;
	//AKS = ref->AKS;
	AKL = ref->AKL;
	DCH = ref->DCH;
	HOD = ref->HOD;
	HAC = ref->HAC;
	LKR = ref->LKR;
	NHO = ref->NHO;
	MUV = ref->MUV;
	MBX = ref->MBX;
	NTR = ref->NTR;
	LV3 = ref->LV3;
	LV3Trig = ref->LV3Trig;
	LV3TrigRare = ref->LV3TrigRare;
	LV3ABTrig = ref->LV3ABTrig;
	LV3ATrigRare = ref->LV3ATrigRare;
	LV3BTrigRare = ref->LV3BTrigRare;
	for(int i=0; i<10; i++) ChTrEff[i] = ref->ChTrEff[i];

	return *this;
}
