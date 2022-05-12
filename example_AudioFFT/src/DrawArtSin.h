/************************************************************
************************************************************/
#pragma once

/************************************************************
************************************************************/
#include <ofMain.h>
#include <ofxAudioFFT.h>
#include <ofxVboSet.h>
#include "sj_common.h"

/************************************************************
************************************************************/

/**************************************************
**************************************************/
class DRAW_ARTSIN{
private:
	/****************************************
	****************************************/
	ofx__VBO_SET VboSet_N;
	ofx__VBO_SET VboSet_A;
	ofx__VBO_SET VboSet_MixDown;
	
	ofTrueTypeFont font_L;
	ofTrueTypeFont font_S;
	
	/****************************************
	****************************************/
	void Refresh_vboVerts(const AUDIO_FFT& fft, bool b_Overlay = false);
	void Refresh_vboColor();
	
public:
	DRAW_ARTSIN();
	~DRAW_ARTSIN();
	
	void setup(const AUDIO_FFT& fft);
	void update(const AUDIO_FFT& fft);
	void draw(ofFbo& fbo);
};

