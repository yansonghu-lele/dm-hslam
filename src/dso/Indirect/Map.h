#pragma once

#include <boost/thread.hpp>
#include "util/settings.h"
#include "util/NumType.h"

#include <memory>
#include <set>

namespace dso
{
    class MapPoint;
    class Frame;

    class Map
    {
    private:
        std::set<std::shared_ptr<MapPoint>> mspMapPoints;
        std::set<std::shared_ptr<Frame>> mspKeyFrames;

        std::vector<std::shared_ptr<MapPoint>> mvpReferenceMapPoints;

        size_t mnMaxMPid;
        size_t mnMaxKFid;

        // Index related to a big change in the map (loop closure, global BA)
        int mnBigChangeIdx;

        boost::mutex mMutexMap;

        bool Busy = false; // is pose graph running?
        boost::mutex mutexBusy;

    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        Map();

        ~Map();
        void AddKeyFrame(std::shared_ptr<Frame> pKF);
        void AddMapPoint(std::shared_ptr<MapPoint> pMP);
        void EraseMapPoint(std::shared_ptr<MapPoint> pMP);
        void EraseKeyFrame(std::shared_ptr<Frame> pKF);
        void SetReferenceMapPoints(const std::vector<std::shared_ptr<MapPoint>> &vpMPs);
        void InformNewBigChange();
        int GetLastBigChangeIdx();

        void GetAllKeyFrames(std::vector<std::shared_ptr<Frame>> &_Out);
        void GetAllMapPoints(std::vector<std::shared_ptr<MapPoint>> &_Out);
        std::vector<std::shared_ptr<MapPoint>> GetReferenceMapPoints();

        long unsigned int MapPointsInMap();
        long unsigned KeyFramesInMap();

        size_t GetMaxMPid();

        void clear();

        bool isIdle()
        {
            boost::unique_lock<boost::mutex> lock(mutexBusy);
            return !Busy;
        }

        void setBusy(bool _isBusy)
        {
            boost::unique_lock<boost::mutex> lock(mutexBusy);
            Busy = _isBusy;
        }

        std::vector<std::shared_ptr<Frame>> mvpKeyFrameOrigins;
        boost::mutex mMutexMapUpdate;
    };
} // namespace dso