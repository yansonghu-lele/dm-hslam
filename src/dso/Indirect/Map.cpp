#include "Indirect/Map.h"

#include "Indirect/MapPoint.h"
#include "Indirect/Frame.h"
#include "util/globalFuncs.h"
#include "util/FrameShell.h"

namespace dso
{
    using namespace std;

    Map::Map()
    {
        mnMaxKFid = 0;
        mnMaxMPid = 0;
        mnBigChangeIdx = 0;
    }

    Map::~Map()
    {
        clear();
    }

    void Map::AddKeyFrame(shared_ptr<Frame> pKF)
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mspKeyFrames.insert(pKF);
        if (pKF->fs->KfId > mnMaxKFid)
            mnMaxKFid = pKF->fs->KfId;
    }

    void Map::AddMapPoint(shared_ptr<MapPoint> pMP)
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mspMapPoints.insert(pMP);
        if ( pMP->id > mnMaxMPid)
            mnMaxMPid = pMP->id;
    }

    void Map::EraseMapPoint(shared_ptr<MapPoint> pMP)
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mspMapPoints.erase(pMP);
    }

    void Map::EraseKeyFrame(shared_ptr<Frame> pKF)
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mspKeyFrames.erase(pKF);
    }

    void Map::SetReferenceMapPoints(const vector<shared_ptr<MapPoint>> &vpMPs)
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mvpReferenceMapPoints = vpMPs;
    }

    void Map::InformNewBigChange()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        mnBigChangeIdx++;
    }

    int Map::GetLastBigChangeIdx()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        return mnBigChangeIdx;
    }

    void Map::GetAllKeyFrames(vector<shared_ptr<Frame>>& _Out)
    {
        releaseVec(_Out);
        boost::lock_guard<boost::mutex> l(mMutexMap);
        _Out.reserve(mspKeyFrames.size());
        _Out.assign(mspKeyFrames.begin(), mspKeyFrames.end());
    }

    void Map::GetAllMapPoints(vector<shared_ptr<MapPoint>>& _Out)
    {
        releaseVec(_Out);
        boost::lock_guard<boost::mutex> l(mMutexMap);
        _Out.reserve(mspMapPoints.size());
        _Out.assign(mspMapPoints.begin(), mspMapPoints.end());
    }

    long unsigned int Map::MapPointsInMap()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        return mspMapPoints.size();
    }

    long unsigned int Map::KeyFramesInMap()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        return mspKeyFrames.size();
    }

    vector<shared_ptr<MapPoint>> Map::GetReferenceMapPoints()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        return mvpReferenceMapPoints;
    }

    size_t Map::GetMaxMPid()
    {
        boost::lock_guard<boost::mutex> l(mMutexMap);
        return mnMaxMPid;
    }


    void Map::clear()
    {
        mspMapPoints.clear();

        mspMapPoints.clear();
        mspKeyFrames.clear();
        mnMaxKFid = 0;
        mvpReferenceMapPoints.clear();
        mvpKeyFrameOrigins.clear();
    }

    
} // namespace dso