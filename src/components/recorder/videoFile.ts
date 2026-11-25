/**
 * The file a MediaRecorder writes is a live stream: a fragmented MP4 (or a WebM) whose header may give no duration,
 * at a variable frame rate (a tab capture sends no frame while nothing moves). Players cope; upload sites often do
 * not: LinkedIn, for one, reads the header and refuses the video as shorter than 3 s. The recording is turned here,
 * in the browser, into a plain MP4, the index (moov) first with its duration:
 * - re-encoded (WebCodecs) to H.264 at a constant 30 frames per second where the browser can;
 * - else remuxed: the same frames copied into a plain MP4, which needs no codec at all.
 */
import { MP4, WEBM, BlobSource, BufferTarget, Conversion, Input, Mp4OutputFormat, Output, QUALITY_HIGH, canEncodeVideo, type ConversionVideoOptions } from 'mediabunny';

/** Largest frame kept: within LinkedIn's limit (4096 × 2304) and what H.264 level 5.1 encodes everywhere. */
const MAX_W = 3840, MAX_H = 2160;
const FPS = 30;

export interface PlainVideo { blob: Blob; how: 'encoded' | 'remuxed' }

async function convert(blob: Blob, video: (w: number, h: number) => Promise<ConversionVideoOptions | null>, onProgress?: (f: number) => void): Promise<Blob | null> {
  const input = new Input({ source: new BlobSource(blob), formats: [MP4, WEBM] });
  try {
    const track = await input.getPrimaryVideoTrack();
    if (!track) return null;
    const opts = await video(track.displayWidth, track.displayHeight);
    if (!opts) return null;
    const output = new Output({ format: new Mp4OutputFormat({ fastStart: 'in-memory' }), target: new BufferTarget() });
    const conversion = await Conversion.init({ input, output, showWarnings: false, video: opts, audio: { discard: true } });
    if (!conversion.isValid) return null;
    if (onProgress) conversion.onProgress = f => onProgress(f);
    await conversion.execute();
    const buf = (output.target as BufferTarget).buffer;
    return buf ? new Blob([buf], { type: 'video/mp4' }) : null;
  } finally {
    input.dispose();
  }
}

/** H.264 at a constant frame rate, the frame scaled to fit, its sides even (H.264 codes 2 × 2 chroma blocks). */
const encoded = async (w: number, h: number): Promise<ConversionVideoOptions | null> => {
  if (typeof VideoEncoder === 'undefined') return null;
  const k = Math.min(1, MAX_W / w, MAX_H / h);
  const even = (v: number) => Math.max(2, 2 * Math.floor(v * k / 2));
  const width = even(w), height = even(h);
  if (!(await canEncodeVideo('avc', { width, height }))) return null;
  return { codec: 'avc', width, height, fit: 'contain', frameRate: FPS, quality: QUALITY_HIGH, keyFrameInterval: 2 };
};

/**
 * A plain MP4 of the recording `blob`, or null where it can be neither re-encoded nor remuxed (it is then kept as
 * recorded). `onProgress` gets a fraction between 0 and 1.
 */
export async function toPlainMp4(blob: Blob, onProgress?: (f: number) => void): Promise<PlainVideo | null> {
  try {
    const out = await convert(blob, encoded, onProgress);
    if (out) return { blob: out, how: 'encoded' };
  } catch { /* an encoder that fails midway: remuxed below */ }
  const out = await convert(blob, async () => ({}), onProgress);
  return out ? { blob: out, how: 'remuxed' } : null;
}
