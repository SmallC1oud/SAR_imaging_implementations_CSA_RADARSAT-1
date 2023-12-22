# SAR_imaging_implementations
This project focuses on Synthetic Aperture Radar (SAR) imaging algorithms and implementations. In the first stage, we will realize four parts: point target simulation and matched filtering, CS imaging algorithm, RD imaging algorithm, and implement a motion compensation algorithm; I will add other contents if I have time.

本次实验使用的是RADARSAT-1原始数据，该数据采集于2002年6月16日，照射的是加拿大温哥华地区。数据储存在 dat_01.001 文件中，包含近19400条记录，每八条记录包含一条传输脉冲的复制信号。每条距离线有9288个复回波采样点，按照uint类型存储。除去复制信号，该记录共18818字节，先是192字节头信息和50字节辅助信息，然后是18576字节回波数据。

# 使用CSA对星载SAR数据进行成像，下图是成像实例

（其斜距分辨率为6m，地面分辨率为10m，单视对全方位带宽进行处理，相应的分辨率为9m）
![图片](https://github.com/SmallC1oud/SAR_imaging_implementations/assets/77475570/03ab34f4-135e-4195-967e-29f2caccf363)


注意：需要先读取RADARSAT-1数据，然后修改路径，由于数据过大，这里只上传了CDdata1.mat作为测试，完整数据可访问 https://www.alipan.com/s/61jtFLwLz6T 进行下载

# 读取数据说明
1、打开specify_parameters.m文件

2、修改输入输出路径（文件中的路径是我电脑的路径，需要自行修改）

3、运行specify_parameters.m文件得到两个SAR回波数据（可修改其中的参数得到不同区域的数据，进行成像）

4、数据在CD_run_params.mat中有储存
详细说明见Extracting_SAR_Data_from_the_CD.pdf文件

# 结果展示
原始数据：
![图片](https://github.com/SmallC1oud/SAR_imaging_implementations/assets/77475570/4888357c-3954-4f57-a347-1b02eb6938cf)

经RC、SRC、一致RCMC后，在距离多普勒域数据已被拉直：
![图片](https://github.com/SmallC1oud/SAR_imaging_implementations/assets/77475570/5387c0df-defe-4cdb-ad56-65714cb71b8d)

CSA成像结果：
![图片](https://github.com/SmallC1oud/SAR_imaging_implementations/assets/77475570/d45d832e-5c03-4a03-be4c-f386614c080e)

成像灰度图：
![图片](https://github.com/SmallC1oud/SAR_imaging_implementations/assets/77475570/67046f60-11c1-42ce-8c3d-44556f2d31fd)
