function [Y, lbls, Ytest, lblstest] = loadLocalData(dataset, baseDir, dirSep)
% LOADLOCALDATA 
% SHEFFIELDML

lbls = [];
lblstest = [];
Ytest=[];


if length(dataset) > 6 && strcmp(dataset(1:7), 'gesture')
    load([baseDir 'gestures/mat/' dataset]);
    Y=data;
    lbls = [height ;width];
    
    return
end

try   
    switch dataset
        %--------- LOCALLY STORED--
        case 'weizmann_horses'
            [Ytr, height, width] = datasetFromImages([localDatasetsDirectoryLarge 'weizmann_horses/']);
            load([baseDir 'weizmann_horses/weizmann_horses_test.mat']);
            Y = Ytr;
            lbls = [height, width];
            Ytest = Yts;
        case 'mulan_emotions'
            load([baseDir 'mulan/emotions/mulan_emotions.mat']);
            lblstest = lblsTest;
        case 'humanPoseAll'
            load([baseDir 'humanPose/humanPoseAll.mat']);
            Yall{1} = Y;
            Yall{2} = Yim;
            Yall{3} = Z;
            YtestAll{1} = Y_test;
            YtestAll{2} = Yim_test;
            YtestAll{3} = Z_test;
            Y = Yall;
            Ytest = YtestAll;
            clear lbls
            lbls = {height,width, seq};
        case 'YaleFace04'
            load([baseDir 'Yale/YaleFace04.mat']);
            lbls = [height;width];
        case 'YaleSubset6_1'
            load([baseDir 'Yale/YaleSubset6.mat']);
            lbls = [height;width];
        case 'YaleSubset4_2'
            load([baseDir 'Yale/YaleSubset4_2.mat']);
            lbls = [height;width];
        case 'YaleSubset4_1'
            load([baseDir 'Yale/subset4.mat']);
            Y{1} = Y02; Y{2}=Y04; Y{3}=Y26; Y{4}=Y37;
            lbls = [192;168];
        case 'humanPose'
            load([baseDir 'sharedVargplvm/humanPose.mat']);
            Ytmp = Y; clear 'Y';
            Y{1}=Ytmp;Y{2}=Y_test;Y{3}=Z;Y{4}=Z_test;
        case 'eth-car1'
            %load('data/car1_masked');
            load([baseDir 'sharedVargplvm/car1_masked']);
            lbls = [256; 256];
        case 'eth-horse5'
            load([baseDir 'sharedVargplvm/horse5_masked']);
            lbls = [256; 256];
        case 'eth-tomato1'
            load([baseDir 'sharedVargplvm/tomato1_masked']);
            lbls = [256; 256];
        case 'ocean'
            load([baseDir 'video/DATA_Ocean.mat']);
        case 'missa'
            load([baseDir 'video' dirSep 'miss-americaHD.mat']);
            lbls = [288 ;360];
        case 'horse'
            load([baseDir 'video/DATA_Horse.mat']);
            lbls = [height ;width];
        case 'horse2'
            load([baseDir 'video/DATA_Horse2.mat']);
            lbls = [height ;width];
        case 'horse2cropped'
            load([baseDir 'video/DATA_Horse2cropped.mat']);
            lbls = [height ;width];
        case 'dog'
            load([baseDir 'video/DATA_DogTr']);
            lbls = [height ;width];
        case 'dogFilled'
            load([baseDir 'video/DATA_DogTrFilled']);
            lbls = [height ;width];
        case 'dog2'
            load([baseDir 'video/DATA_DogTr2']);
            lbls = [height ;width];
        case 'dino'
            load([baseDir 'video/DATA_dino']);
            lbls = [height ;width];
        case 'dinoFilled'
            load([baseDir 'video/DATA_dinoFilled']);
            lbls = [height ;width];
        case 'claire'
            load([baseDir 'video/DATA_claire']);
            lbls = [height ;width];
        case 'susie'
            load([baseDir 'video/DATA_susie']);
            lbls = [height ;width];
        case 'grandma'
            load([baseDir 'video/DATA_grandma']);
            Y = Y(700:800,:);
            lbls = [height ;width];
        case 'ADN'
            load([baseDir 'video/DATA_ADN']);
            lbls = [height ;width];
        case 'head'
            load([baseDir 'video/DATA_head']);
            Y = Y(15:135,:);
            lbls = [height ;width];
        case 'USecon'
            load([baseDir 'finance/DATA_USecon']);
        case 'NYSE2Small'
            load([baseDir 'finance/nyse_new2Small']);
            lbls = dataSetInfo;
        otherwise
           % load([baseDir dirSep dataset]);
            load([baseDir dataset]);
    end
catch e
    throw(e);
end