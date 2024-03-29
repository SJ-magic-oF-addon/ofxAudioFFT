/************************************************************
■参考
	of_v0.11.2_osx_release/examples/sound/audioInputExample
	of_v0.11.2_osx_release/examples/sound/audioOutputExample


■Audio I/F
input:cable
	モノラル。
	
	walkmanからinputする時は、LR(白赤)に分けて入れる。
	ちなみに、「標準ステレオ」をinputしたつもりで動作させると、
	「低音スカスカなんだけど」って感じになった。
	
input:Level
	walkmanから入れる場合、cable上のLevelはLine Level。
	input感度は、minimum〜設定していく。
	
■primitive study
	■DTM初心者のためのオーディオインターフェイス選び 2016
		http://www.dtmstation.com/archives/51971048.html
	
	■オーディオ信号レベル
		http://detail.chiebukuro.yahoo.co.jp/qa/question_detail/q1413957214
			マイク 1.6mV
			エレキギター 200mV
			ライン 1,230mV（1,000mV）
			ヘッドフォン、スピーカー 2,400-2.830mV
			
	■４つのオーディオ信号レベル【音楽・動画などの制作・配信に必須】
		https://masafumiiwasaki.com/blog/audio-signal-level-4types/
		
	■バランスとアンバランス【オーディオ機器のケーブル接続】
		https://masafumiiwasaki.com/blog/balance-unbalance/#toc6
	

■【今さら聞けない用語シリーズ】
	■初心者向け【今さら聞けない用語シリーズ】初めてのオーディオ・インターフェース選び
		http://info.shimamura.co.jp/digital/knowledge/2014/06/21591
		
	■【今さら聞けない用語シリーズ】3分でわかる！インピーダンスって何でしょう？ ロー出しハイ受け？
		http://info.shimamura.co.jp/digital/knowledge/2014/01/16805
n		
	■【今さら聞けない用語シリーズ】3分でわかる！オーディオ・ケーブルの種類とバランス・アンバランスの違い
		http://info.shimamura.co.jp/digital/knowledge/2014/01/16501
		
	■【今さら聞けない用語シリーズ】音と周波数
		http://info.shimamura.co.jp/digital/knowledge/2014/03/19260
		
	■【今さら聞けない用語シリーズ】3分でわかる！コンデンサーマイクとダイナミックマイクの違い おすすめマイク
		http://info.shimamura.co.jp/digital/knowledge/2014/04/23017
************************************************************/
#pragma once

/************************************************************
************************************************************/
#include <ofMain.h>
#include <ofxAudioFFT.h>

#include "sj_common.h"

#include "DrawFFT.h"
#include "DrawArtSin.h"

/************************************************************
************************************************************/

class ofApp : public ofBaseApp{
private:
	/****************************************
	****************************************/
	int AUDIO_IN_CHS = 2;
	int AUDIO_OUT_CHS = 2;
	
	ofSoundStream *soundStream;
	bool b_EnableAudioOut;
	int soundStream_Input_DeviceId;
	int soundStream_Output_DeviceId;
	
	ofFbo fbo_out;
	ofImage img_Frame;
	
	AUDIO_FFT AudioFFT;
	
	/********************
	********************/
	DRAW_FFT DrawFFT;
	DRAW_ARTSIN DrawArtSin;
	
	/****************************************
	****************************************/
	void Sound__Start();
	void Sound__Stop();
	void Sound__Reset();
	void setup_Gui();
	void set_CamParam();
	void draw_frame(ofFbo& fbo);
	void update_AudioFFT_ParamForFFTCal();
	
public:
	/****************************************
	****************************************/
	
	/****************************************
	****************************************/
	ofApp(int _soundStream_Input_DeviceId, int _soundStream_Output_DeviceId, int _AUDIO_IN_CHS, int _AUDIO_OUT_CHS);
	~ofApp();
	
	void setup_sound();
	
	void setup() override;
	void update() override;
	void draw() override;
	
	void exit() override;
	
	void audioIn(ofSoundBuffer & buffer) override;
	void audioOut(ofSoundBuffer & buffer) override;
	
	void keyPressed(int key) override;
	void keyReleased(int key) override;
	void mouseMoved(int x, int y ) override;
	void mouseDragged(int x, int y, int button) override;
	void mousePressed(int x, int y, int button) override;
	void mouseReleased(int x, int y, int button) override;
	void mouseEntered(int x, int y) override;
	void mouseExited(int x, int y) override;
	void windowResized(int w, int h) override;
	void dragEvent(ofDragInfo dragInfo) override;
	void gotMessage(ofMessage msg) override;
	
};
