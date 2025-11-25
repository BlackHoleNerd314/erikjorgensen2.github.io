import numpy as np
import pygame
from TTS.api import TTS


def float32_to_int16(audio):
    # Clip to [-1,1] just in case, then scale float32 to int16 range
    audio = np.clip(audio, -1, 1)
    return (audio * 32767).astype(np.int16)


def main():
    pygame.mixer.init(frequency=11025, size=-16, channels=1)  # match TTS sample rate and mono

    tts = TTS(model_name="tts_models/en/ljspeech/tacotron2-DDC", gpu=False)

    text = "Hello! This is Coqui TTS playing audio directly through pygame without saving."

    # Generate waveform numpy array
    wav = tts.tts(text)

    # Convert float32 to int16 PCM
    pcm_data = float32_to_int16(wav)

    # Convert to bytes
    pcm_bytes = pcm_data.tobytes()

    # Create pygame Sound object from raw data
    sound = pygame.mixer.Sound(buffer=pcm_bytes)

    # Play it
    sound.play()

    # Keep program running until playback ends
    while pygame.mixer.get_busy():
        pygame.time.wait(100)


if __name__ == "__main__":
    main()
