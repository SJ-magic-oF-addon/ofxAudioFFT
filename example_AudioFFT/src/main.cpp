#include "ofMain.h"
#include "ofApp.h"

//========================================================================
int main( int argc, char** argv ){
	
	/********************
	********************/
	//Use ofGLFWWindowSettings for more options like multi-monitor fullscreen
	ofGLWindowSettings settings;
	settings.setSize(1024, 768);
	settings.windowMode = OF_WINDOW; //can also be OF_FULLSCREEN
	
	auto window = ofCreateWindow(settings);
	
	/********************
	********************/
	int AUDIO_IN_CHS = 2;
	int AUDIO_OUT_CHS = 2;
	
	int soundStream_Input_DeviceId = -1;
	int soundStream_Output_DeviceId = -1;
	
	/********************
	********************/
	printf("---------------------------------\n");
	printf("> parameters\n");
	printf("\t-i      Audio in (-1)\n");
	printf("\t-o      Audio out(-1)\n");
	printf("\t-ichs   Audio in  chs(2)\n");
	printf("\t-ochs   Audio out chs(2)\n");
	printf("---------------------------------\n");
	
	for(int i = 1; i < argc; i++){
		if( strcmp(argv[i], "-i") == 0 ){
			if(i+1 < argc) { soundStream_Input_DeviceId = atoi(argv[i+1]); i++; }
		}else if( strcmp(argv[i], "-o") == 0 ){
			if(i+1 < argc) { soundStream_Output_DeviceId = atoi(argv[i+1]); i++; }
		}else if( strcmp(argv[i], "-ichs") == 0 ){
			if(i+1 < argc) { AUDIO_IN_CHS = atoi(argv[i+1]); i++; }
		}else if( strcmp(argv[i], "-ochs") == 0 ){
			if(i+1 < argc) { AUDIO_OUT_CHS = atoi(argv[i+1]); i++; }
		}
	}
	
	/********************
	********************/
	printf("> parameters\n");
	printf("sound_In  = %d\n", soundStream_Input_DeviceId);
	printf("sound_Out = %d\n", soundStream_Output_DeviceId);
	fflush(stdout);
	
	/********************
	********************/
	ofRunApp(window, make_shared<ofApp>( soundStream_Input_DeviceId, soundStream_Output_DeviceId, AUDIO_IN_CHS, AUDIO_OUT_CHS ));
	ofRunMainLoop();
}


