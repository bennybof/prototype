/* main.js
   Expects:
   - data/rules.json (in GitHub site)
   - Audio stems + secret_versions hosted at rules.assets.audio_base_url
   - UI elements in HTML with IDs:
       #btnDownloadWav, #btnDownloadMp3, #btnRefresh
       #statusText (optional)
       #visualLayer (optional container for popups)
   - Optional: background video element #bgVideo (optional)
   - Optional: include lamejs in HTML for MP3 export:
       <script src="https://cdn.jsdelivr.net/npm/lamejs@1.2.1/lame.min.js"></script>
*/

(() => {
  // -------------------------
  // Utility: status + logging
  // -------------------------
  const $ = (sel) => document.querySelector(sel);
  const statusEl = $("#statusText");
  function setStatus(msg) {
    console.log(msg);
    if (statusEl) statusEl.textContent = msg;
  }

  // -------------------------
  // Time parsing: "mm:ss:cs" -> seconds
  // -------------------------
  function parseTimeMMSSCS(t) {
    // Accept "MM:SS:CS" (centiseconds), also allow "HH:MM:SS:CS" or "SS:CS" loosely if needed
    if (typeof t !== "string") throw new Error("Time must be string, got: " + typeof t);
    const parts = t.split(":").map(p => p.trim());
    if (parts.length < 2) throw new Error("Bad time format: " + t);

    let cs = 0, s = 0, m = 0, h = 0;

    if (parts.length === 2) {
      // "SS:CS"
      s = parseInt(parts[0], 10);
      cs = parseInt(parts[1], 10);
    } else if (parts.length === 3) {
      // "MM:SS:CS"
      m = parseInt(parts[0], 10);
      s = parseInt(parts[1], 10);
      cs = parseInt(parts[2], 10);
    } else if (parts.length === 4) {
      // "HH:MM:SS:CS"
      h = parseInt(parts[0], 10);
      m = parseInt(parts[1], 10);
      s = parseInt(parts[2], 10);
      cs = parseInt(parts[3], 10);
    } else {
      throw new Error("Bad time format: " + t);
    }

    if ([h, m, s, cs].some(n => Number.isNaN(n))) throw new Error("Bad time numeric: " + t);
    return h * 3600 + m * 60 + s + (cs / 100);
  }

  // -------------------------
  // Seeded PRNG: Mulberry32
  // -------------------------
  function hashStringToUint32(str) {
    // FNV-1a style
    let h = 2166136261 >>> 0;
    for (let i = 0; i < str.length; i++) {
      h ^= str.charCodeAt(i);
      h = Math.imul(h, 16777619);
    }
    return h >>> 0;
  }

  function mulberry32(seed) {
    let a = seed >>> 0;
    return function rand() {
      a |= 0;
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function randBetween(rand, min, max) {
    return min + (max - min) * rand();
  }

  function chance(rand, p) {
    return rand() < p;
  }

  // -------------------------
  // Fetch + decode audio
  // -------------------------
  async function fetchArrayBuffer(url) {
  const res = await fetch(url);
  const ct = res.headers.get("content-type") || "(no content-type)";
  console.log("[fetch]", res.status, ct, url);

  if (!res.ok) {
    const text = await res.text().catch(() => "");
    throw new Error(`Fetch failed ${res.status} for ${url} (content-type: ${ct}) ${text.slice(0, 80)}`);
  }

  return await res.arrayBuffer();
}

async function decodeToBuffer(audioCtxLike, url) {
  try {
    const arr = await fetchArrayBuffer(url);
    return await audioCtxLike.decodeAudioData(arr.slice(0));
  } catch (e) {
    throw new Error(`DECODE FAILED for ${url} :: ${e.message || e}`);
  }
}

  // -------------------------
  // WAV encoding (PCM 16-bit)
  // -------------------------
  function audioBufferToWavBlob(buffer, bitDepth = 16) {
    // Only 16-bit PCM implemented here (simple + widely compatible)
    if (bitDepth !== 16) throw new Error("Only 16-bit WAV supported in this encoder.");
    const numCh = buffer.numberOfChannels;
    const sampleRate = buffer.sampleRate;
    const numFrames = buffer.length;

    // interleave
    const interleaved = new Float32Array(numFrames * numCh);
    for (let ch = 0; ch < numCh; ch++) {
      const data = buffer.getChannelData(ch);
      for (let i = 0; i < numFrames; i++) {
        interleaved[i * numCh + ch] = data[i];
      }
    }

    // PCM16
    const bytesPerSample = 2;
    const blockAlign = numCh * bytesPerSample;
    const byteRate = sampleRate * blockAlign;
    const dataSize = numFrames * blockAlign;

    const bufferSize = 44 + dataSize;
    const view = new DataView(new ArrayBuffer(bufferSize));

    // RIFF header
    writeString(view, 0, "RIFF");
    view.setUint32(4, 36 + dataSize, true);
    writeString(view, 8, "WAVE");

    // fmt chunk
    writeString(view, 12, "fmt ");
    view.setUint32(16, 16, true);        // PCM chunk size
    view.setUint16(20, 1, true);         // format = PCM
    view.setUint16(22, numCh, true);
    view.setUint32(24, sampleRate, true);
    view.setUint32(28, byteRate, true);
    view.setUint16(32, blockAlign, true);
    view.setUint16(34, 16, true);        // bits per sample

    // data chunk
    writeString(view, 36, "data");
    view.setUint32(40, dataSize, true);

    // samples
    let offset = 44;
    for (let i = 0; i < interleaved.length; i++, offset += 2) {
      let s = interleaved[i];
      // clamp
      s = Math.max(-1, Math.min(1, s));
      // convert
      view.setInt16(offset, s < 0 ? s * 0x8000 : s * 0x7FFF, true);
    }

    return new Blob([view.buffer], { type: "audio/wav" });

    function writeString(dv, start, str) {
      for (let i = 0; i < str.length; i++) dv.setUint8(start + i, str.charCodeAt(i));
    }
  }

  // -------------------------
  // MP3 encoding (optional, requires lamejs)
  // -------------------------
  function audioBufferToMp3Blob(buffer, kbps = 320) {
    if (typeof lamejs === "undefined") {
      throw new Error("MP3 export requires lamejs. Add it to HTML via a <script> tag.");
    }
    const numCh = buffer.numberOfChannels;
    const sampleRate = buffer.sampleRate;

    // lamejs wants Int16 samples
    const mp3Encoder = new lamejs.Mp3Encoder(numCh, sampleRate, kbps);

    const blockSize = 1152;
    const mp3Data = [];

    const chData = [];
    for (let ch = 0; ch < numCh; ch++) chData.push(buffer.getChannelData(ch));

    const numSamples = buffer.length;

    function floatTo16BitPCM(float32) {
      const out = new Int16Array(float32.length);
      for (let i = 0; i < float32.length; i++) {
        let s = Math.max(-1, Math.min(1, float32[i]));
        out[i] = s < 0 ? (s * 0x8000) : (s * 0x7FFF);
      }
      return out;
    }

    // encode in blocks
    for (let i = 0; i < numSamples; i += blockSize) {
      const end = Math.min(i + blockSize, numSamples);

      if (numCh === 1) {
        const mono = floatTo16BitPCM(chData[0].subarray(i, end));
        const mp3buf = mp3Encoder.encodeBuffer(mono);
        if (mp3buf.length) mp3Data.push(new Uint8Array(mp3buf));
      } else {
        const left = floatTo16BitPCM(chData[0].subarray(i, end));
        const right = floatTo16BitPCM(chData[1].subarray(i, end));
        const mp3buf = mp3Encoder.encodeBuffer(left, right);
        if (mp3buf.length) mp3Data.push(new Uint8Array(mp3buf));
      }
    }

    const flush = mp3Encoder.flush();
    if (flush.length) mp3Data.push(new Uint8Array(flush));

    return new Blob(mp3Data, { type: "audio/mpeg" });
  }

  // -------------------------
  // Download helper
  // -------------------------
  function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 5000);
  }

  // -------------------------
  // Visual layer helper (very simple)
  // -------------------------
  const visualLayer = $("#visualLayer");
  function popupImage({ file, holdSeconds = 2 }) {
    if (!visualLayer) return;
    const img = document.createElement("img");
    img.src = file;
    img.style.position = "absolute";
    img.style.left = "50%";
    img.style.top = "50%";
    img.style.transform = "translate(-50%, -50%)";
    img.style.maxWidth = "70vw";
    img.style.maxHeight = "70vh";
    img.style.opacity = "0";
    img.style.transition = "opacity 250ms ease";
    visualLayer.appendChild(img);
    requestAnimationFrame(() => (img.style.opacity = "1"));
    setTimeout(() => {
      img.style.opacity = "0";
      setTimeout(() => img.remove(), 400);
    }, Math.max(200, holdSeconds * 1000));
  }

  // -------------------------
  // UI random placement
  // -------------------------
  function rectsOverlap(a, b) {
    return !(a.right <= b.left || a.left >= b.right || a.bottom <= b.top || a.top >= b.bottom);
  }

  function randomPlaceButtons(rules, rand) {
    const wavBtn = $("#btnDownloadWav");
    const mp3Btn = $("#btnDownloadMp3");
    const refreshBtn = $("#btnRefresh");
    const btns = [wavBtn, mp3Btn, refreshBtn].filter(Boolean);
    if (!btns.length) return;

    const rp = rules?.ui?.buttons?.random_placement;
    if (!rp?.enabled) return;

    const margin = rp.margin_px ?? 24;
    const avoidTopFrac = rp.avoid_top_fraction ?? 0.15;
    const maxAttempts = rp.max_attempts ?? 50;

    // Ensure buttons are positioned absolute (CSS should do this too)
    for (const b of btns) {
      b.style.position = "absolute";
    }

    const vw = window.innerWidth;
    const vh = window.innerHeight;

    // Place one by one with collision-avoid
    const placedRects = [];

    for (const b of btns) {
      // temporarily set visibility to measure
      b.style.visibility = "hidden";
      b.style.left = "0px";
      b.style.top = "0px";
      b.style.visibility = "visible";

      const br = b.getBoundingClientRect();
      const bw = br.width;
      const bh = br.height;

      const minX = margin;
      const maxX = Math.max(margin, vw - bw - margin);
      const minY = Math.max(margin + vh * avoidTopFrac, margin);
      const maxY = Math.max(minY, vh - bh - margin);

      let placed = false;
      for (let attempt = 0; attempt < maxAttempts; attempt++) {
        const x = Math.floor(randBetween(rand, minX, maxX));
        const y = Math.floor(randBetween(rand, minY, maxY));
        const rect = { left: x, top: y, right: x + bw, bottom: y + bh };

        const overlaps = placedRects.some(r => rectsOverlap(rect, r));
        if (!overlaps) {
          b.style.left = `${x}px`;
          b.style.top = `${y}px`;
          placedRects.push(rect);
          placed = true;
          break;
        }
      }

      if (!placed) {
        // fallback: stack near bottom
        const idx = placedRects.length;
        const x = margin;
        const y = Math.floor(vh - (idx + 1) * (bh + 12) - margin);
        b.style.left = `${x}px`;
        b.style.top = `${Math.max(margin, y)}px`;
        placedRects.push({ left: x, top: y, right: x + bw, bottom: y + bh });
      }
    }
  }

  // -------------------------
  // Rules loading
  // -------------------------
  async function loadRules() {
    const res = await fetch("data/rules.json", { cache: "no-store" });
    if (!res.ok) throw new Error("Failed to load data/rules.json");
    return await res.json();
  }

  // -------------------------
  // Core: build render plan from rules
  // -------------------------
  function decideSeed() {
    // New seed each refresh: use crypto if available for uniqueness
    if (crypto?.getRandomValues) {
      const buf = new Uint32Array(1);
      crypto.getRandomValues(buf);
      return buf[0] >>> 0;
    }
    return (Date.now() >>> 0);
  }

  function makeRandFromSeed(seed) {
    return mulberry32(seed >>> 0);
  }

  function buildSelectionPlan(rules, seed) {
    const rand = makeRandFromSeed(seed);

    // 1) Secret bypass?
    const secrets = rules?.secrets?.bypass_versions ?? [];
    let chosenSecret = null;
    for (const s of secrets) {
      const p = s?.probability ?? 0;
      if (p > 0 && chance(rand, p)) {
        chosenSecret = s;
        break;
      }
    }

    // 2) Intro selection
    const introRule = rules?.structure?.intro_rule;
    const introSelected = introRule ? chance(rand, introRule.roll?.probability ?? 0) : false;

    // 3) FX flags
    const fx = rules?.fx ?? {};
    const convolution = fx?.convolution_reverb_vox;
    const convolutionEnabled = convolution ? chance(rand, convolution.roll?.probability ?? 0) : false;

    // 4) Stem inclusion decisions
    const stems = rules?.stems ?? [];

    const stemGlobalIncluded = new Map();   // stemId -> bool
    const stemPerLoopIncluded = new Map();  // stemId -> bool[] length = full_loops_count (+ possible extra)
    const perLoopState = new Map();         // for checking "selected in loop"

    const fullLoopsCount = rules?.structure?.full_loops_count ?? 0;

    // First pass: global selection (scope=global)
    for (const st of stems) {
      const id = st.id;
      // bundles
      const hasFilesArray = Array.isArray(st.files);
      const roll = st.roll || {};
      let include = true;

      if (roll.scope === "global") {
        include = chance(rand, roll.probability ?? 0);
      } else if (roll.scope === "per_full_loop") {
        // not globally gated unless you want it
        include = true;
      } else {
        include = true;
      }

      // requires (global)
      if (Array.isArray(st.requires) && st.requires.length) {
        // temporarily mark; we’ll enforce after we compute globals
      }

      stemGlobalIncluded.set(id, include);

      // pre-init per loop arrays
      stemPerLoopIncluded.set(id, Array(fullLoopsCount).fill(false));
      perLoopState.set(id, Array(fullLoopsCount).fill(false));

      // bundles: if global false, none
      if (hasFilesArray) {
        // keep only by id; per trigger scheduling uses files
      }
    }

    // Enforce requires: if required stem not globally included, force false
    for (const st of stems) {
      const id = st.id;
      const req = st.requires || [];
      if (req.length) {
        const ok = req.every(rid => stemGlobalIncluded.get(rid) === true);
        if (!ok) stemGlobalIncluded.set(id, false);
      }
    }

    // Second pass: per-loop selection
    for (const st of stems) {
      const id = st.id;
      const roll = st.roll || {};
      const globalGate = stemGlobalIncluded.get(id) === true;

      const perLoop = stemPerLoopIncluded.get(id);

      if (!globalGate) {
        // never active
        continue;
      }

      if (roll.scope === "per_full_loop") {
        for (let li = 0; li < fullLoopsCount; li++) {
          perLoop[li] = chance(rand, roll.probability ?? 0);
        }
      } else {
        // default: if global selected, treat as active for all loops unless there is if_selected_then_roll
        for (let li = 0; li < fullLoopsCount; li++) perLoop[li] = true;
      }

      // if_selected_then_roll: overrides per loop activeness
      if (st.if_selected_then_roll && globalGate) {
        const p = st.if_selected_then_roll.probability ?? 0;
        const scope = st.if_selected_then_roll.scope;
        if (scope === "per_full_loop") {
          for (let li = 0; li < fullLoopsCount; li++) {
            perLoop[li] = chance(rand, p);
          }
        }
      }
    }

    // Extra segment rule (adds 1 extra loop segment, and shifts vox_section_2)
    const extraRule = rules?.structure?.extra_segment_rule;
    let extraInserted = false;
    if (extraRule?.roll?.probability) {
      extraInserted = chance(rand, extraRule.roll.probability);
    }

    // Build final plan object
    return {
      seed,
      rand, // keep for UI placement, etc.
      chosenSecret,
      introSelected,
      convolutionEnabled,
      extraInserted,
      stemGlobalIncluded,
      stemPerLoopIncluded
    };
  }

  // -------------------------
  // Schedule events based on triggers
  // -------------------------
  function computeTimeline(rules, plan) {
    const fullLoopDur = parseTimeMMSSCS(rules.structure.full_loop_duration);
    const miniLoopDur = parseTimeMMSSCS(rules.structure.mini_loop_duration);
    const introEnd = parseTimeMMSSCS(rules.structure.intro_end_time);

    const loopsStart = plan.introSelected ? introEnd : 0;
    const fullLoopsCount = rules.structure.full_loops_count;

    let totalDuration = loopsStart + fullLoopsCount * fullLoopDur;

    // extra segment inserted after 4.5 loops into loopsStart timeline
    let extraStart = null;
    if (plan.extraInserted) {
      extraStart = loopsStart + (4 * fullLoopDur) + (0.5 * fullLoopDur);
      totalDuration += fullLoopDur; // extend total duration by one full loop
    }

    return { fullLoopDur, miniLoopDur, introEnd, loopsStart, totalDuration, extraStart };
  }

  // -------------------------
  // LFO curve generation (simple automation curve)
  // -------------------------
  function makeSmoothedRandomLFO(rand, rateHz, durationSec, samplePointsPerSec = 100) {
    // step at rateHz but smooth between steps
    const points = Math.max(2, Math.floor(durationSec * samplePointsPerSec));
    const out = new Float32Array(points);

    const stepInterval = 1 / rateHz;           // seconds per step
    let nextStepTime = 0;
    let currentTarget = randBetween(rand, -1, 1);
    let prevTarget = currentTarget;

    for (let i = 0; i < points; i++) {
      const t = (i / (points - 1)) * durationSec;

      if (t >= nextStepTime) {
        prevTarget = currentTarget;
        currentTarget = randBetween(rand, -1, 1);
        nextStepTime += stepInterval;
      }

      // Smooth within current step window
      const windowStart = nextStepTime - stepInterval;
      const alpha = stepInterval > 0 ? Math.max(0, Math.min(1, (t - windowStart) / stepInterval)) : 1;
      // Smoothstep
      const sm = alpha * alpha * (3 - 2 * alpha);
      out[i] = prevTarget + (currentTarget - prevTarget) * sm;
    }

    return out;
  }

  function makeSineLFO(rateHz, durationSec, samplePointsPerSec = 100) {
    const points = Math.max(2, Math.floor(durationSec * samplePointsPerSec));
    const out = new Float32Array(points);
    for (let i = 0; i < points; i++) {
      const t = (i / (points - 1)) * durationSec;
      out[i] = Math.sin(2 * Math.PI * rateHz * t);
    }
    return out;
  }

  function averageCurves(a, b) {
    const n = Math.min(a.length, b.length);
    const out = new Float32Array(n);
    for (let i = 0; i < n; i++) out[i] = (a[i] + b[i]) / 2;
    return out;
  }

  // -------------------------
  // Rendering
  // -------------------------
  async function renderVersion(rules, plan, outputFormat) {
    // If secret bypass: just download the file (or transcode if mp3 requested)
    if (plan.chosenSecret) {
      const secretUrl = rules.assets.audio_base_url.replace(/\/$/, "") + "/" + plan.chosenSecret.file;
      setStatus(`Secret version selected. Fetching ${outputFormat.toUpperCase()}...`);
      const arr = await fetchArrayBuffer(secretUrl);

      if (outputFormat === "wav") {
        return { blob: new Blob([arr], { type: "audio/wav" }), filename: `secret_${plan.seed}.wav` };
      }

      // MP3 requested
      // Decode WAV -> encode MP3 (needs lamejs)
      const tempCtx = new (window.AudioContext || window.webkitAudioContext)();
      const buf = await tempCtx.decodeAudioData(arr.slice(0));
      await tempCtx.close?.();
      const mp3Blob = audioBufferToMp3Blob(buf, rules.timebase.mp3_bitrate_kbps ?? 320);
      return { blob: mp3Blob, filename: `secret_${plan.seed}.mp3` };
    }

    setStatus("Loading audio + building timeline...");
    const timeline = computeTimeline(rules, plan);

    // Offline context
    const sampleRate = rules.timebase.sample_rate_hz ?? 44100;
    const lengthFrames = Math.ceil(timeline.totalDuration * sampleRate);
    const offline = new OfflineAudioContext(2, lengthFrames, sampleRate);

    // Master bus
    const masterGain = offline.createGain();
    masterGain.gain.value = 1.0;
    masterGain.connect(offline.destination);

    // Intro fade-in if intro not selected
    if (!plan.introSelected) {
      const fade = rules?.structure?.loops_start_rule?.if_intro_not_selected_master_fade_in;
      if (fade?.type === "linear_gain_ramp") {
        const dur = parseTimeMMSSCS(fade.duration);
        masterGain.gain.setValueAtTime(fade.from_gain ?? 0, 0);
        masterGain.gain.linearRampToValueAtTime(fade.to_gain ?? 1, Math.max(0.001, dur));
      }
    }

    // Optional convolution FX bus
    let convolver = null, convWet = null;
    if (plan.convolutionEnabled) {
      const fx = rules.fx.convolution_reverb_vox;
      const irUrl = rules.assets.audio_base_url.replace(/\/$/, "") + "/" + fx.impulse_response;
      setStatus("Loading impulse response...");
      const irBuf = await decodeToBuffer(offline, irUrl);

      convolver = offline.createConvolver();
      convolver.buffer = irBuf;

      convWet = offline.createGain();
      // wet_gain_db -> linear
      const wetDb = fx.wet_gain_db ?? -20;
      convWet.gain.value = Math.pow(10, wetDb / 20);

      convolver.connect(convWet);
      convWet.connect(masterGain);
    }

    // LFO (simple automation curves)
    const lfoSpec = rules?.fx?.lfo;
    const lfoByStemParam = new Map(); // key `${stemId}:${param}` -> curve Float32Array

    if (lfoSpec?.sources?.length && Array.isArray(lfoSpec.targets)) {
      const duration = timeline.totalDuration;
      const sine = makeSineLFO(lfoSpec.sources[0].rate_hz ?? 0.25, duration, 100);
      const smrd = makeSmoothedRandomLFO(plan.rand, lfoSpec.sources[1].rate_hz ?? 0.25, duration, 100);
      let curve = averageCurves(sine, smrd);

      const [minOut, maxOut] = lfoSpec.output_range ?? [-0.2, 0.2];
      // scale curve from [-1,1] -> [minOut,maxOut]
      const mid = (minOut + maxOut) / 2;
      const amp = (maxOut - minOut) / 2;
      for (let i = 0; i < curve.length; i++) curve[i] = mid + amp * curve[i];

      for (const tgt of lfoSpec.targets) {
        lfoByStemParam.set(`${tgt.stem_id}:${tgt.param}`, curve);
      }
    }

    // Load + schedule stems
    const stems = rules.stems ?? [];
    const audioBase = rules.assets.audio_base_url.replace(/\/$/, "");

    // For conditional rules like pre_snare gain depending on oh_rev selection per loop,
    // we need to know per-loop selections of oh_rev.
    const perLoopIncluded = plan.stemPerLoopIncluded;

    setStatus("Loading stem buffers...");
    const bufferCache = new Map(); // url -> AudioBuffer

    async function getBufByPath(path) {
      const url = audioBase + "/" + path;
      if (bufferCache.has(url)) return bufferCache.get(url);
      const buf = await decodeToBuffer(offline, url);
      bufferCache.set(url, buf);
      return buf;
    }

    function scheduleSource(buf, whenSec, stemId, filePath) {
      const src = offline.createBufferSource();
      src.buffer = buf;

      // per-stem chain
      const g = offline.createGain();
      g.gain.value = 1.0;

      // optional pan node (for LFO / fixed pan)
      const panner = offline.createStereoPanner();
      panner.pan.value = 0;

      src.connect(g);
      g.connect(panner);

      // Apply convolution only for vox stems when enabled
      const shouldConv = plan.convolutionEnabled && (filePath.includes("vox") || filePath.includes("Vox") || filePath.includes("VOX"));
      if (shouldConv && convolver) {
        // split dry/wet
        const dry = offline.createGain();
        dry.gain.value = 1.0;
        panner.connect(dry);
        dry.connect(masterGain);

        const wetSend = offline.createGain();
        wetSend.gain.value = 1.0;
        panner.connect(wetSend);
        wetSend.connect(convolver);
      } else {
        panner.connect(masterGain);
      }

      // Apply LFO automation if specified
      const gainCurve = lfoByStemParam.get(`${stemId}:gain`);
      if (gainCurve) {
        // gain is additive mapping: base 1 + curve
        // We'll apply to gain param using setValueCurveAtTime around 1 baseline:
        // create curve for automation: 1 + curve[i]
        const auto = new Float32Array(gainCurve.length);
        for (let i = 0; i < gainCurve.length; i++) auto[i] = 1 + gainCurve[i];
        g.gain.setValueCurveAtTime(auto, 0, timeline.totalDuration);
      }
      const panCurve = lfoByStemParam.get(`${stemId}:pan`);
      if (panCurve) {
        // additive mapping: base 0 + curve
        panner.pan.setValueCurveAtTime(panCurve, 0, timeline.totalDuration);
      }

      src.start(Math.max(0, whenSec));
    }

    // Apply special pan flip for vlins targets
    if (rules?.fx?.global_pan_flip_for_vlins) {
      const flip = rules.fx.global_pan_flip_for_vlins;
      const flipOn = chance(plan.rand, flip.roll?.probability ?? 0);
      const panVal = flipOn ? (flip.if_true_pan ?? 0.64) : (flip.else_pan ?? -0.64);
      // We'll store a fixed pan automation for these stems by inserting a constant curve
      for (const targetId of (flip.targets ?? [])) {
        const curve = new Float32Array(2);
        curve[0] = panVal;
        curve[1] = panVal;
        lfoByStemParam.set(`${targetId}:pan`, curve);
      }
    }

    // Visual triggers during decision time (NOT during render)
    tryTriggerDecisionVisuals(rules, plan);

    // Schedule intro if selected
    if (plan.introSelected && rules?.structure?.intro_rule?.file) {
      const buf = await getBufByPath(rules.structure.intro_rule.file);
      scheduleSource(buf, 0, rules.structure.intro_rule.id ?? "intro_section", rules.structure.intro_rule.file);
    }

    // Extra segment: shift vox_section_2 if inserted
    const vox2Shift = (plan.extraInserted ? parseTimeMMSSCS("00:17:16") : 0);

    // Helper: compute event times for triggers
    function eventTimeFromTrigger(tr, loopIndex) {
      const loopStart = timeline.loopsStart + loopIndex * timeline.fullLoopDur;
      switch (tr.type) {
        case "full_loop_start":
          return loopStart;
        case "full_loop_offset":
          return loopStart + parseTimeMMSSCS(tr.offset);
        case "mini_loop_start": {
          // mini_loops are 1-based
          // We schedule at each mini loop listed
          return loopStart + (0 * timeline.miniLoopDur); // caller handles list
        }
        default:
          return loopStart;
      }
    }

    // After_full_loops type triggers (global timeline after loopsStart)
    function timeAfterFullLoops(fullLoops, offsetStr) {
      return timeline.loopsStart + (fullLoops * timeline.fullLoopDur) + (offsetStr ? parseTimeMMSSCS(offsetStr) : 0);
    }

    // Schedule standard stems
    setStatus("Scheduling stems...");
    const fullLoopsCount = rules.structure.full_loops_count;

    for (const st of stems) {
      const id = st.id;

      // Skip if global not included
      if (plan.stemGlobalIncluded.get(id) !== true) continue;

      const perLoop = plan.stemPerLoopIncluded.get(id) || Array(fullLoopsCount).fill(false);

      // Bundle case: st.files array
      if (Array.isArray(st.files)) {
        // all_or_none bundles: if included, schedule each file using same triggers
        const files = st.files;
        for (let li = 0; li < fullLoopsCount; li++) {
          if (!perLoop[li]) continue;
          for (const tr of (st.triggers ?? [])) {
            if (tr.type === "full_loop_start") {
              const t = timeline.loopsStart + li * timeline.fullLoopDur;
              for (const f of files) {
                const buf = await getBufByPath(f);
                scheduleSource(buf, t, id, f);
              }
            }
          }
        }
        continue;
      }

      // Standard single file stem
      const file = st.file;

      // Special global triggers not tied to per-loop triggers
      const triggers = st.triggers ?? [];

      for (const tr of triggers) {
        if (tr.type === "after_full_loops") {
          const t = timeAfterFullLoops(tr.full_loops ?? 0, tr.offset);
          const buf = await getBufByPath(file);
          // apply shift if vox_section_2 and extra inserted
          const shiftedT = (id === "vox_section_2") ? (t + vox2Shift) : t;
          scheduleSource(buf, shiftedT, id, file);
          continue;
        }

        if (tr.type === "after_full_loops_plus") {
          // Use what rules.json encoded
          const base = timeAfterFullLoops(tr.full_loops ?? 0, null);
          const plusMini = (tr.plus_mini_loops ?? 0) * timeline.miniLoopDur;
          const off = tr.offset ? parseTimeMMSSCS(tr.offset) : 0;
          const t = base + plusMini + off;
          const buf = await getBufByPath(file);
          scheduleSource(buf, t, id, file);
          continue;
        }

        // Per-loop triggers
        for (let li = 0; li < fullLoopsCount; li++) {
          if (!perLoop[li]) continue;

          // Conditional modifier example: pre_snare gain if oh_rev NOT selected in loop
          // Implemented by pre-adjusting gain via a simple gain multiplier (3 dB)
          // We apply by creating a separate file path? We'll instead apply at scheduling time:
          // (We do it below by checking st.conditional_modifiers)
          const loopStart = timeline.loopsStart + li * timeline.fullLoopDur;

          if (tr.type === "mini_loop_start") {
            const minis = tr.mini_loops ?? [];
            for (const ml of minis) {
              const t = loopStart + (Math.max(1, ml) - 1) * timeline.miniLoopDur;
              const buf = await getBufByPath(file);

              // schedule with optional gain modifier
              scheduleWithConditionalGain(buf, t, st, li, id, file);
            }
            continue;
          }

          if (tr.type === "full_loop_start") {
            const t = loopStart;
            const buf = await getBufByPath(file);
            scheduleWithConditionalGain(buf, t, st, li, id, file);
            continue;
          }

          if (tr.type === "full_loop_offset") {
            const t = loopStart + parseTimeMMSSCS(tr.offset);
            const buf = await getBufByPath(file);
            scheduleWithConditionalGain(buf, t, st, li, id, file);
            continue;
          }
        }
      }
    }

    // Extra segment events from extra_rule.when_selected.events
    if (plan.extraInserted && timeline.extraStart != null) {
      const extraEvents = rules?.structure?.extra_segment_rule?.when_selected?.events ?? [];
      for (const ev of extraEvents) {
        if (ev.type === "play_stem_once") {
          const t = timeline.extraStart + parseTimeMMSSCS(ev.offset ?? "00:00:00");
          const buf = await getBufByPath(ev.file);
          scheduleSource(buf, t, "secret_lyric", ev.file);
        }
      }
    }

    // Render
    setStatus("Rendering (offline)...");
    const rendered = await offline.startRendering();

    // Encode
    if (outputFormat === "wav") {
      setStatus("Encoding WAV...");
      const wavBlob = audioBufferToWavBlob(rendered, rules.timebase.wav_bit_depth ?? 16);
      return { blob: wavBlob, filename: `render_${plan.seed}.wav` };
    }

    setStatus("Encoding MP3...");
    const mp3Blob = audioBufferToMp3Blob(rendered, rules.timebase.mp3_bitrate_kbps ?? 320);
    return { blob: mp3Blob, filename: `render_${plan.seed}.mp3` };

    // ---- local helper ----
    function scheduleWithConditionalGain(buf, whenSec, st, loopIndex, stemId, filePath) {
      // Base schedule but allow +3dB when condition met
      // We implement by wrapping scheduleSource with gain bump if needed.
      // Simple approach: duplicate scheduleSource logic here with custom gain.
      const src = offline.createBufferSource();
      src.buffer = buf;

      const g = offline.createGain();
      let gainMul = 1.0;

      const mods = st.conditional_modifiers ?? [];
      for (const m of mods) {
        if (m?.apply?.type === "gain_db" && m?.when?.stem_not_selected_in_loop) {
          const depId = m.when.stem_not_selected_in_loop;
          const depLoop = perLoopIncluded.get(depId) || [];
          const depSelected = depLoop[loopIndex] === true;
          if (!depSelected) {
            const db = m.apply.value_db ?? 0;
            gainMul *= Math.pow(10, db / 20);
          }
        }
      }
      g.gain.value = gainMul;

      const panner = offline.createStereoPanner();
      panner.pan.value = 0;

      src.connect(g);
      g.connect(panner);

      const shouldConv = plan.convolutionEnabled && (filePath.includes("vox") || filePath.includes("Vox") || filePath.includes("VOX"));
      if (shouldConv && convolver) {
        const dry = offline.createGain();
        dry.gain.value = 1.0;
        panner.connect(dry);
        dry.connect(masterGain);

        const wetSend = offline.createGain();
        wetSend.gain.value = 1.0;
        panner.connect(wetSend);
        wetSend.connect(convolver);
      } else {
        panner.connect(masterGain);
      }

      // LFO automation if present
      const gainCurve = lfoByStemParam.get(`${stemId}:gain`);
      if (gainCurve) {
        const auto = new Float32Array(gainCurve.length);
        for (let i = 0; i < gainCurve.length; i++) auto[i] = gainMul * (1 + gainCurve[i]);
        g.gain.setValueCurveAtTime(auto, 0, timeline.totalDuration);
      }
      const panCurve = lfoByStemParam.get(`${stemId}:pan`);
      if (panCurve) {
        panner.pan.setValueCurveAtTime(panCurve, 0, timeline.totalDuration);
      }

      src.start(Math.max(0, whenSec));
    }
  }

  function tryTriggerDecisionVisuals(rules, plan) {
    const vts = rules?.ui?.visual_triggers ?? [];
    for (const vt of vts) {
      if (vt?.when?.fx_flag_enabled === "convolution_reverb_vox") {
        if (plan.convolutionEnabled) {
          const action = vt.action;
          if (action?.type === "popup_image") {
            popupImage({ file: action.file, holdSeconds: action.hold_duration_seconds ?? 2 });
          }
        }
      }
      if (vt?.when?.stem_selected === "secret_lyric") {
        // only if extraInserted (secret_lyric is played)
        if (plan.extraInserted) {
          // If you have a background element, you could swap it here
          // For now: show a quick popup
          popupImage({ file: vt.action?.file, holdSeconds: 2 });
        }
      }
    }
  }

  // -------------------------
  // UI SFX (real-time) for button clicks
  // -------------------------
  let realtimeCtx = null;
  async function playUiSfx(rules, plan) {
    // Optional: play UI sounds from GitHub-hosted site paths in rules.ui.on_download_pressed
    const sfxRules = rules?.ui?.on_download_pressed?.play_sounds ?? [];
    if (!sfxRules.length) return;

    if (!realtimeCtx) realtimeCtx = new (window.AudioContext || window.webkitAudioContext)();

    // We do lightweight decoding in realtime context
    for (const entry of sfxRules) {
      if (entry.file) {
        // direct play
        const url = entry.file; // from GitHub site
        if (chance(plan.rand, entry.probability ?? 1.0)) {
          const buf = await decodeToBuffer(realtimeCtx, url);
          const src = realtimeCtx.createBufferSource();
          src.buffer = buf;
          const g = realtimeCtx.createGain();
          const db = entry.gain_db ?? 0;
          g.gain.value = Math.pow(10, db / 20);
          src.connect(g);
          g.connect(realtimeCtx.destination);
          src.start();
        }
      } else if (entry.choose_one && Array.isArray(entry.choose_one)) {
        // choose one by probability weights
        let r = plan.rand();
        let cum = 0;
        for (const opt of entry.choose_one) {
          cum += (opt.probability ?? 0);
          if (r <= cum) {
            const buf = await decodeToBuffer(realtimeCtx, opt.file);
            const src = realtimeCtx.createBufferSource();
            src.buffer = buf;
            src.connect(realtimeCtx.destination);
            src.start();
            break;
          }
        }
      }
    }
  }

  // -------------------------
  // App init
  // -------------------------
  let rules = null;
  let currentSeed = null;
  let currentPlan = null;

  async function refreshPlan() {
    if (!rules) rules = await loadRules();

    currentSeed = decideSeed();
    currentPlan = buildSelectionPlan(rules, currentSeed);

    // Place buttons deterministically by seed
    randomPlaceButtons(rules, makeRandFromSeed(currentSeed));

    setStatus(`Ready. Seed: ${currentSeed}${currentPlan.chosenSecret ? " (SECRET!)" : ""}`);
  }

  async function onDownload(format) {
    try {
      if (!rules) rules = await loadRules();
      if (!currentPlan) await refreshPlan();

      // UI click sound + decision-time visuals
      await playUiSfx(rules, currentPlan);

      const { blob, filename } = await renderVersion(rules, currentPlan, format);
      downloadBlob(blob, filename);

      setStatus(`Downloaded ${filename}`);
    } catch (err) {
      console.error(err);
      setStatus(`Error: ${err.message || err}`);
      alert(err.message || String(err));
    }
  }

  // Wire up UI
  async function init() {
    try {
      rules = await loadRules();
      await refreshPlan();

      const wavBtn = $("#btnDownloadWav");
      const mp3Btn = $("#btnDownloadMp3");
      const refreshBtn = $("#btnRefresh");

      if (wavBtn) wavBtn.addEventListener("click", () => onDownload("wav"));
      if (mp3Btn) mp3Btn.addEventListener("click", () => onDownload("mp3"));
      if (refreshBtn) refreshBtn.addEventListener("click", refreshPlan);

      // If MP3 encoder not present, inform user
      if (typeof lamejs === "undefined" && mp3Btn) {
        mp3Btn.title = "MP3 export requires lamejs. Add the lamejs script tag in index.html.";
      }

      setStatus("Ready.");
    } catch (err) {
      console.error(err);
      setStatus(`Init error: ${err.message || err}`);
    }
  }

  // Start
  window.addEventListener("load", init);

})();

