#include <vector>
#include <TVector3.h>

#ifndef __EXPORTCLASSES__
#define __EXPORTCLASSES__

struct vtxtracks;
struct superTimeOffset;
struct superBurst;
struct superCmpEvent;
struct SCvertex;
struct trak;
struct cluster;

//##### COMPLETE
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

/*class NAbcog_params{
public:
	NAbcog_params():
		alpha(0), alpha_coeff(0), beta(0), beta_coeff(0), mkp(0), mkperr(0),
		mkn(0), mknerr(0), cogX1p(0), cogY1p(0), cogX1n(0), cogY1n(0),
		cogX4p(0), cogY4p(0), cogX4n(0), cogY4n(0), status(0), pkp(0),
		pkdxdzp(0), pkdydzp(0), pkxoffp(0), pkyoffp(0), pkm(0),
		pkdxdzm(0), pkdydzm(0), pkxoffm(0), pkyoffm(0){};

	NAbcog_params(abcog_params_t &ref);
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
};
*/

class NCluster : public TObject{
public:
	NCluster():time(0), dDeadCell(0), E(0){};
	NCluster(cluster &ref);
public:
	float time;
	float dDeadCell;
	/*float x;
	float y;*/
	float E;
	TVector3 position;

	ClassDefNV(NCluster, 1);
};

class NTrak : public TObject{
public:
	NTrak():q(999), vtxID(-1), time(0), p(0), bdxdz(0), bdydz(0){};
	NTrak(trak &ref);
public:
	int q;
	int vtxID;
	float time;
	float p;
	/*float bx;
	float by;*/
	float bdxdz;
	float bdydz;
	TVector3 bDetPos;
	TVector3 aDetPos;
	TVector3 bMomentum;
	TVector3 aMomentum;

	ClassDefNV(NTrak, 1);
};

class NPhysicsTrack : public TObject{
public:
	NPhysicsTrack():
		trackID(-1), clusterID(-1), p(0), E(0){};
	~NPhysicsTrack(){};
public:
	int trackID;
	int clusterID;
	float p;
	float E;
	TVector3 momentum;

	ClassDefNV(NPhysicsTrack, 1);
};

class NPhysicsCluster : public TObject{
public:
	int clusterID;
	float E;
	TVector3 position;

	ClassDefNV(NPhysicsCluster, 1);
};

class NSCVertex : public TObject{
public:
	NSCVertex():Nvtxtrack(0), charge(999), cda(0), chi2(-1){};
	NSCVertex(SCvertex &ref);
public:
	unsigned int Nvtxtrack;
	int charge;
	float cda;
	/*float x;
	float y;
	float z;*/
	float chi2;
	TVector3 position;
	std::vector<NVtxTrack> vtxtrack;

	ClassDefNV(NSCVertex, 1);
};

//class NDCHTrack{
//	float pq;
//	float p;
//	float q;
//	float perr;
//	float chi2;
//	float bx;
//	float by;
//	float bdxdz;
//	float bdydz;
//	float x;
//	float y;
//	float dxdz;
//	float dydz;
//	float time;
//	float quality;
//	float hodtime;
//	int hodstatus;
//	int nhits;
//	unsigned int Nhit;
//	float exhac;
//	float eyhac;
//	float ddeadcell;
//	float sigxx;
//	float sigyy;
//
//	float sigdxdx;
//	float sigdydy;
//	float sigxdx;
//	float sigxy;
//	float sigdxy;
//	float sigxdy;
//	float sigdxdy;
//	float sigydy;
//	int HitPattern;
//	int efficiency[2];
//	int spareInt[2];
//	int spareFloat[2];
//	int LKRclu;
//	float EovP;
//	float Espy;
//	float muvTime;
//	float anavar[20];
//	int anaflag[5];
//};

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

	ClassDefNV(ROOTRawEvent, 1);
};

class ROOTCorrectedEvent : public TObject{
public:
	ROOTCorrectedEvent():
		weight(0){};

	void clear(){
		weight = 0;

		pTrack.clear();
		pCluster.clear();
		goodTracks.clear();
	};

public:
	float weight;

	std::vector<NPhysicsTrack> pTrack;
	std::vector<NPhysicsCluster> pCluster;

	std::vector<int> goodTracks;

	ClassDefNV(ROOTCorrectedEvent, 1);
};

/*class PhysicsEvent : public TObject{

	ClassDefNV(PhysicsEvent, 1);
};*/

class ROOTBurst : public TObject{
public:
	ROOTBurst():nrun(-1), time(-1){};

	ROOTBurst& operator=(superBurst *ref);

	void clear(){
		nrun = -1;
		time = -1;
	};
public:
	int nrun;
	int time;
	NSuperTimeOffset tOffst;

	ClassDefNV(ROOTBurst, 1);
};

#endif
