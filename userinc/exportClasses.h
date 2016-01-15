#ifndef __EXPORTCLASSES__
#define __EXPORTCLASSES__

#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

struct vtxtracks;
struct superTimeOffset;
struct superBurst;
struct superCmpEvent;
struct SCvertex;
struct trak;
struct cluster;
struct GeomCompact;
struct DETstatus;

//##########################
//###   RawObjects
//##########################
class NVtxTrack : public TObject{
public:
	NVtxTrack(): iTrack(-1),bdxdz(0),bdydz(0){};
	NVtxTrack(vtxtracks &ref);
	~NVtxTrack(){};

public:
	int iTrack;
	float bdxdz;
	float bdydz;

	ClassDefNV(NVtxTrack, 1);
};

class NCluster : public TObject{
public:
	NCluster():time(0), dDeadCell(0), E(0), lkr_acc(false), iTrack(-99){};
	NCluster(cluster &ref);
public:
	float time;
	float dDeadCell;
	float E;
	bool lkr_acc;
	int iTrack;
	TVector3 position;

	ClassDefNV(NCluster, 1);
};

class NTrak : public TObject{
public:
	NTrak():p(0), q(999), quality(0), chi2(0), by(0), bx(0), bdxdz(0), bdydz(0), vtxID(-1), time(0), dDeadCell(0){};
	NTrak(trak &ref);
public:
	float p;
	int q;
	float quality;
	float chi2;
	float by;
	float bx;
	float bdxdz;
	float bdydz;
	int vtxID;
	float time;
	float dDeadCell;
	TVector3 bDetPos;
	TVector3 aDetPos;
	TVector3 bMomentum;
	TVector3 aMomentum;

	ClassDefNV(NTrak, 1);
};

//##########################
//###   Corrected objects
//##########################
class NPhysicsTrack : public TObject{
public:
	NPhysicsTrack():
		trackID(-1), clusterID(-1), p(0), E(0), lkr_acc(false){};
	~NPhysicsTrack(){};
public:
	int trackID;
	int clusterID;
	float p;
	float E;
	bool lkr_acc;
	TVector3 momentum;

	ClassDefNV(NPhysicsTrack, 1);
};

class NPhysicsCluster : public TObject{
public:
	NPhysicsCluster(): clusterID(-1), E(0){};
	~NPhysicsCluster(){};
public:
	int clusterID;
	float E;
	TVector3 position;

	ClassDefNV(NPhysicsCluster, 1);
};

class NSCVertex : public TObject{
public:
	NSCVertex():Nvtxtrack(0), charge(999), cda(0), chi2(-1), time(-1){};
	NSCVertex(SCvertex &ref);
public:
	unsigned int Nvtxtrack;
	int charge;
	float cda;
	float chi2;
	float time;
	TVector3 position;
	std::vector<NVtxTrack> vtxtrack;

	ClassDefNV(NSCVertex, 1);
};

//##########################
//###   Physics Objects
//##########################
class NRecoParticle : public TObject{
public:
	NRecoParticle(): pdgID(0), parentTrack(-1), parentCluster(-1), parentVertex(-1){};
	NRecoParticle(int id): pdgID(id), parentTrack(-1), parentCluster(-1), parentVertex(-1){};
	~NRecoParticle(){};

	void clear(){
		parentTrack = -1;
		parentCluster = -1;
		parentVertex = -1;
		vertex.SetXYZ(0,0,0);
		P.SetXYZM(0,0,0,0);
	};
public:
	int pdgID;
	int parentTrack;
	int parentCluster;
	int parentVertex;
	TVector3 vertex;
	TLorentzVector P;

	ClassDefNV(NRecoParticle, 1);
};

class NMCParticle : public TObject{
public:
	NMCParticle(): pdgID(0){};
	NMCParticle(int id): pdgID(id){};
	~NMCParticle(){};

	void clear(){
		vertex.SetXYZ(0,0,0);
		P.SetXYZM(0,0,0,0);
	};
public:
	int pdgID;
	TVector3 vertex;
	TLorentzVector P;

	ClassDefNV(NMCParticle, 1);
};


//##########################
//###   Databases
//##########################
class NSuperTimeOffset : public TObject{
public:
	NSuperTimeOffset():
		Version(0), Tag(0), Aks(0),	Kab(0),	Nmv(0),	Akl(0),
		Hod(0),	Nho(0),	Dch(0),	Lkr(0),	Hac(0),	Muv(0),
		LkrTag(0), LkrNhod(0), LkrAkl(0), LkrHac(0),
		KabPlus(0),	KabMinus(0){};

	NSuperTimeOffset& operator=(superTimeOffset &ref);
	~NSuperTimeOffset(){};

public:
	float Version;
	float Tag;
	float Aks;
	float Kab;
	float Nmv;
	float Akl;
	float Hod;
	float Nho;
	float Dch;
	float Lkr;
	float Hac;
	float Muv;
	float LkrTag;
	float LkrNhod;
	float LkrAkl;
	float LkrHac;
	float KabPlus;
	float KabMinus;

	ClassDefNV(NSuperTimeOffset, 1);
};

class NAbcog_params{
public:
	NAbcog_params():
		alpha(0), alpha_coeff(0), beta(0), beta_coeff(0), mkp(0), mkperr(0),
		mkn(0), mknerr(0), cogX1p(0), cogY1p(0), cogX1n(0), cogY1n(0),
		cogX4p(0), cogY4p(0), cogX4n(0), cogY4n(0), status(0), pkp(0),
		pkdxdzp(0), pkdydzp(0), pkxoffp(0), pkyoffp(0), pkm(0),
		pkdxdzm(0), pkdydzm(0), pkxoffm(0), pkyoffm(0){};

	NAbcog_params& operator=(void *ref);
public:
	float alpha;
	float alpha_coeff;
	float beta;
	float beta_coeff;
	float mkp;
	float mkperr;
	float mkn;
	float mknerr;
	float cogX1p;
	float cogY1p;
	float cogX1n;
	float cogY1n;
	float cogX4p;
	float cogY4p;
	float cogX4n;
	float cogY4n;
	int   status;
	float pkp;
	float pkdxdzp;
	float pkdydzp;
	float pkxoffp;
	float pkyoffp;
	float pkm;
	float pkdxdzm;
	float pkdydzm;
	float pkxoffm;
	float pkyoffm;

	ClassDefNV(NAbcog_params, 1);
};

class Nxyz : public TObject{
public:
	Nxyz(): x(0), y(0), z(0){};
	~Nxyz(){};

	void SetXYZ(double x, double y, double z){
		this->x = x;
		this->y = y;
		this->z = z;
	}
public:
	double x, y, z;

	ClassDefNV(Nxyz, 1);
};

class NDCH : public TObject{
public:
	NDCH(){};
	~NDCH(){};
public:
	Nxyz PosChamber;

	ClassDefNV(NDCH, 1);
};

class NGeom : public TObject{
public:
	NGeom(){};
	~ NGeom(){};

	NGeom& operator=(GeomCompact *ref);
public:
	NDCH Dch[4];
	Nxyz Lkr;
	ClassDefNV(NGeom, 1);
};

class NDETStatus : public TObject{
public:
	NDETStatus(){};
	~NDETStatus(){};

	NDETStatus& operator=(DETstatus *ref);
public:
	int TAG;
	int AKS;
	int AKL;
	int DCH;
	int HOD;
	int HAC;
	int LKR;
	int NHO;
	int MUV;
	int MBX;
	int NTR;
	int LV3;
	int LV3Trig;
	int LV3TrigRare;
	int LV3ABTrig;
	int LV3ATrigRare;
	int LV3BTrigRare;
	int ChTrEff[10];
	ClassDefNV(NDETStatus, 1);
};

//##########################
//###   Top nodes
//##########################

class ROOTRawEvent : public TObject{
public:
	ROOTRawEvent():
		Nvtx(0), Ncluster(0), timeStamp(-1), trigWord(0){};

	ROOTRawEvent& operator=(superCmpEvent *ref);

	void clear(){
		Nvtx = 0;
		Ncluster = 0;
		timeStamp = -1;
		trigWord = 0;

		vtx.clear();
		track.clear();
		cluster.clear();
	};

public:
	unsigned int Nvtx;
	unsigned int Ncluster;
	int timeStamp;
	int trigWord;

	std::vector<NSCVertex> vtx;
	std::vector<NTrak> track;
	std::vector<NCluster> cluster;
	NDETStatus DETStatus;

	ClassDefNV(ROOTRawEvent, 1);
};

class ROOTCorrectedEvent : public TObject{
public:
	ROOTCorrectedEvent():
		failedCond(-1), goodVertexID(-1), weight(0), kaonP(0){};

	void clear(){
		weight = 0;
		failedCond = -1;
		goodVertexID = -1;

		kaonP = 0;

		kaonMomentum.SetXYZ(0,0,0);

		pTrack.clear();
		pCluster.clear();
		goodTracks.clear();
	};

public:
	int failedCond;
	int goodVertexID;
	float weight;
	float kaonP;
	TVector3 kaonMomentum;

	std::vector<NPhysicsTrack> pTrack;
	std::vector<NPhysicsCluster> pCluster;

	std::vector<int> goodTracks;

	ClassDefNV(ROOTCorrectedEvent, 1);
};

class ROOTBurst : public TObject{
public:
	ROOTBurst():isData(false), isMC(false), pbWall(false), nrun(-1), time(-1), period(-1), beamCharge(-99), alpha(0){};

	ROOTBurst& operator=(superBurst *ref);

	void clear(){
		isData = false;
		this->isMC = false;
		pbWall = false;
		nrun = -1;
		time = -1;
		period = -1;
		beamCharge = -99;
		alpha = 0;
	};
public:
	bool isData;
	bool isMC;
	bool pbWall;
	int nrun;
	int time;
	int period;
	int beamCharge;
	double alpha;
	NSuperTimeOffset tOffst;
	NAbcog_params abcog_params;

	ClassDefNV(ROOTBurst, 1);
};

class ROOTFileHeader: public TObject{
public:
	ROOTFileHeader():NProcessedEvents(0), NFailedEvents(0), NPassedEvents(0){};
	~ROOTFileHeader(){};

public:
	int NProcessedEvents;
	int NFailedEvents;
	int NPassedEvents;

	ClassDefNV(ROOTFileHeader, 1);
};

class ROOTPhysicsEvent : public TObject{
public:
	ROOTPhysicsEvent(): x(0), y(0), mee(0), pic(), ep(11), em(-11), gamma(22), pi0(111), kaon(), mu(){};
	~ROOTPhysicsEvent(){};

	void clear(){
		x = 0;
		y = 0;
		mee = 0;
		pic.clear();
		ep.clear();
		em.clear();
		gamma.clear();
		pi0.clear();
		kaon.clear();
	};
public:
	double x,y;
	double mee;
	NRecoParticle pic;
	NRecoParticle ep;
	NRecoParticle em;
	NRecoParticle gamma;
	NRecoParticle pi0;
	NRecoParticle kaon;
	NRecoParticle mu;

	ClassDefNV(ROOTPhysicsEvent, 1);
};

class ROOTMCEvent: public TObject{
public:
	ROOTMCEvent(): xTrue(-1), ep(11), em(-11), k(321){};
	~ROOTMCEvent(){};

	void clear(){
		xTrue = -1;
		ep.clear();
		em.clear();
		k.clear();
	}
public:
	float xTrue;
	NMCParticle ep;
	NMCParticle em;
	NMCParticle k;

	ClassDefNV(ROOTMCEvent, 1);
};
#endif
