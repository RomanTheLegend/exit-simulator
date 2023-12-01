#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <thread>
#include <mutex>
#include <chrono>
#include "datapoint.h"
#include "common.h"

#define SLOPE_LEN 4
#define SLEEP_MULTIPLIER 1

class gps
{
public:
    std::string getTimeString() const
    {

        return timestring_;
    }
    long long getTime() const
    {

        return time_;
    }

    double getLat() const
    {

        return lat_;
    }

    double getLon() const
    {

        return lon_;
    }

    double getHMSL() const
    {

        return hMSL_;
    }

    double getVelN() const
    {

        return velN_;
    }

    double getVelE() const
    {

        return velE_;
    }

    double getVelD() const
    {

        return velD_;
    }

    double getHAcc() const
    {

        return hAcc_;
    }

    double getVAcc() const
    {

        return vAcc_;
    }

    double getSAcc() const
    {

        return sAcc_;
    }

    double getHeading() const
    {

        return heading_;
    }

    double getCAcc() const
    {

        return cAcc_;
    }

    double getGPSFix() const
    {

        return gpsFix_;
    }

    double getNumSV() const
    {

        return numSV_;
    }

    std::string timestring_;
    long long time_;
    double lat_;
    double lon_;
    double hMSL_;
    double velN_;
    double velE_;
    double velD_;
    double hAcc_;
    double vAcc_;
    double sAcc_;
    double heading_;
    double cAcc_;
    double gpsFix_;
    double numSV_;
};

class CSVDataReader
{
public:
    CSVDataReader(const std::string &filename) : filename_(filename) {}

    void startReading()
    {
        std::thread readerThread(&CSVDataReader::readData, this);
        readerThread.detach();
    }

    // bool isReady() const
    // {
    //     return isReady_;
    // }

    bool finishReading() const
    {
        return finishReading_;
    }

    // void markAsRead()
    // {
    //     isReady_ = false;
    // }

    gps getGPSdata()
    {
        while (!isReady_)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(2*SLEEP_MULTIPLIER));
        }
        isReady_ = false;
        return dp;
    }

private:
    gps dp = gps();

    std::string filename_;
    bool isReady_ = false;
    bool finishReading_ = false;

    mutable std::mutex mutex_; // For thread safety

    void readData()
    {
        std::ifstream file(filename_);
        if (!file.is_open())
        {
            std::cerr << "Error opening file: " << filename_ << std::endl;
            return;
        }

        std::string line;
        while (std::getline(file, line))
        {
            isReady_ = false;

            std::istringstream iss(line);
            std::string token;

            // Assuming the format of the CSV file is consistent
            std::getline(iss, token, ',');
            dp.timestring_ = token;
            dp.time_ = convertToUnixEpoch(token);

            std::getline(iss, token, ',');
            dp.lat_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.lon_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.hMSL_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.velN_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.velE_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.velD_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.hAcc_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.vAcc_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.sAcc_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.heading_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.cAcc_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.gpsFix_ = std::stod(token);

            std::getline(iss, token, ',');
            dp.numSV_ = std::stod(token);

            isReady_ = true;
            // Notify any potential listeners (not implemented in this example)
            // You can add a callback or observer pattern for notifying other parts of your program
            std::this_thread::sleep_for(std::chrono::milliseconds(20*SLEEP_MULTIPLIER));
        }
        finishReading_ = true;
        file.close();
    }

    long long convertToUnixEpoch(const std::string &timestamp)
    {
        std::tm tm = {};
        std::istringstream iss(timestamp);
        char delimiter;

        // Parse the timestamp string
        iss >> std::get_time(&tm, "%Y-%m-%dT%H:%M:%S%c");

        // Check if there is a fraction of seconds
        if (iss.peek() == '.')
        {
            iss.ignore(); // Ignore the dot separator

            // Read the fraction of seconds
            std::string fractionString;
            std::getline(iss, fractionString, 'Z');

            // Convert fractionString to an integer if it is not empty
            int fractionSeconds = (!fractionString.empty()) ? std::stoi(fractionString) : 0;

            auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));
            auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch());

            // Add the fraction of seconds in milliseconds
            return milliseconds.count() + fractionSeconds * 10 + 3600000;
        }
        else
        {
            auto tp = std::chrono::system_clock::from_time_t(std::mktime(&tm));
            return std::chrono::duration_cast<std::chrono::milliseconds>(tp.time_since_epoch()).count();
        }
    }
};

DataPoint curDp;
DataPoint prevDp;
DataPoint runStartDp;
bool alreadyFalling = false;
long long startTs = 0;
long long exitTs;
DataPoint timeseries[SLOPE_LEN];



double getSlope()
{
    double sumx = 0, sumy = 0, sumxx = 0, sumxy = 0;

    for (int i = 0; i <= SLOPE_LEN-1; ++i)
    {
        double y = timeseries[i].velD;
        double t = timeseries[i].t;

        sumx += t;
        sumy += y;
        sumxx += t * t;
        sumxy += t * y;
    }


    int n = SLOPE_LEN;

    double slope = (sumxy - sumx * sumy / n) / (sumxx - sumx * sumx / n);

    std::cout << " | Slope: " << slope;

    return slope;

}


void addElement(DataPoint value, DataPoint array[], int size)
{
    for (int i = 0; i < size - 1; ++i)
    {
        array[i] = array[i + 1];
    }
    array[size - 1] = value;
}

bool isFreefall()
{
    if (alreadyFalling)
        return true;

    // Get interpolation coefficient
    timeseries[SLOPE_LEN - 1].az = getSlope();
    double g = A_GRAVITY;

int p = SLOPE_LEN -2, c = SLOPE_LEN -1;

    double a = (g - timeseries[p].velD) / (timeseries[c].velD - timeseries[p].velD);
    std::cout << " | a: " << a << std::endl;
    // Check vertical speed
    if (a < 0 || 1 < a)
        return false;

    // Check accuracy
    double vAcc = timeseries[p].vAcc + a * (timeseries[c].vAcc - timeseries[p].vAcc);
    if (vAcc > 10)
        return false;

    // Check acceleration
    double az = timeseries[p].az + a * (timeseries[c].az - timeseries[p].az);
    if (az < g / 5.)
        return false;

    exitTs = (long long)(timeseries[p].ts + a * (timeseries[c].ts - timeseries[p].ts) - g / az*1000.) ;
    return true;
}

int main()
{
    CSVDataReader dataReader("data4.csv");
    dataReader.startReading();

    curDp = DataPoint();
    prevDp = DataPoint();
    int i = 0;

    while (!dataReader.finishReading())
    {
        // Example usage:

        gps gps_ = dataReader.getGPSdata();

        prevDp = curDp;
        curDp.vAcc = gps_.getCAcc();
        curDp.velD = gps_.getVelD();
        curDp.ts = gps_.getTime();

        if (startTs != 0)
        {
            curDp.t = (double)(curDp.ts - startTs) / 1000;
            addElement(curDp, timeseries, SLOPE_LEN);

        }
        else
        {
            startTs = gps_.getTime();
        }

        // std::cout << "Time: " << gps_.getTime() << " | TS: " << curDp.t << " | Latitude: " << gps_.getLat() << " | Longitude: " << gps_.getLon() << " | HMSL: " << gps_.getHMSL() << " | VelN: " << gps_.getVelN() << " | VelE: " << gps_.getVelE() << " | VelD: " << gps_.getVelD() << " | HAcc: " << gps_.getHAcc() << " | VAcc: " << gps_.getVAcc() << " | SAcc: " << gps_.getSAcc() << " | Heading: " << gps_.getHeading() << " | CAcc: " << gps_.getCAcc() << " | GPSFix: " << gps_.getGPSFix() << " | NumSV: " << gps_.getNumSV() << std::endl;
        std::cout << gps_.getTimeString() << " | Timestamp: " << gps_.getTime() << " | Delta t: " << curDp.t  << " | VelD: " << gps_.getVelD() << " | VAcc: " << gps_.getVAcc();

        if (i >= SLOPE_LEN && isFreefall())
        {
            std::cout << "Falling! " << exitTs << std::endl;
            return 0;
        }

        ++i;
    }

    return 0;
}
