import React from "react";
import { AbsoluteFill, Sequence, useCurrentFrame } from "remotion";
import { CellCycle, TwoCells } from "./Cells";
import {
  BLUE, ChartCard, Fonts, Headline, INK, INK2, Kicker, PAPER, Progress, Reveal, Scene, Sub, ORANGE, PURPLE, GREEN,
} from "./theme";

// Paced for reading: every scene holds its text well past the reveal. 1800 frames = 60 s at 30 fps.
const D = {
  hook: 150, gap: 168, mech: 216, memory: 210,
  persist: 210, calib: 225, ident: 210, clinic: 225, close: 186,
};
const order = ["hook", "gap", "mech", "memory", "persist", "calib", "ident", "clinic", "close"] as const;
export const TOTAL = order.reduce((a, k) => a + D[k], 0);

const starts: Record<string, number> = {};
{
  let acc = 0;
  for (const k of order) { starts[k] = acc; acc += D[k]; }
}

const Hook: React.FC = () => (
  <Scene durationInFrames={D.hook}>
    <Reveal><Kicker>the problem</Kicker></Reveal>
    <Reveal delay={4}><Headline>Two cells. Same genome.<br />Same drug.</Headline></Reveal>
    <Reveal delay={16} style={{ marginTop: 30 }}><TwoCells duration={D.hook} /></Reveal>
    <Reveal delay={80}><Sub>One dies. One doesn&rsquo;t. Nothing in the sequence says which.</Sub></Reveal>
  </Scene>
);

const Gap: React.FC = () => (
  <Scene durationInFrames={D.gap}>
    <Reveal><Kicker color={PURPLE}>the gap</Kicker></Reveal>
    <Reveal delay={5}>
      <Headline>Single-cell simulators model<br />regulation &mdash; not growth,<br />division, or inheritance.</Headline>
    </Reveal>
    <Reveal delay={26}>
      <Sub>So they cannot generate the heritable states that decide survival, which is exactly what the methods they benchmark are built to find.</Sub>
    </Reveal>
  </Scene>
);

const Mech: React.FC = () => (
  <Scene durationInFrames={D.mech}>
    <Reveal><Kicker>Biomodelling.jl 2.0</Kicker></Reveal>
    <Reveal delay={4}><Headline>Stochastic kinetics inside<br />growing, dividing cells.</Headline></Reveal>
    <Reveal delay={18} style={{ marginTop: 26 }}><CellCycle duration={D.mech} /></Reveal>
    <Reveal delay={150}>
      <Sub>Genes replicate mid-cycle. Molecules partition at division. Both daughters inherit the promoter state.</Sub>
    </Reveal>
  </Scene>
);

const chartScene = (
  key: keyof typeof D, kicker: string, color: string,
  head: React.ReactNode, src: string, sub: React.ReactNode, maxHeight = 500
) => () => (
  <Scene durationInFrames={D[key]}>
    <Reveal><Kicker color={color}>{kicker}</Kicker></Reveal>
    <Reveal delay={4}><Headline size={60}>{head}</Headline></Reveal>
    <ChartCard src={src} delay={16} maxHeight={maxHeight} />
    <Reveal delay={40}><Sub size={30}>{sub}</Sub></Reveal>
  </Scene>
);

const Memory = chartScene("memory", "what emerges", GREEN,
  <>Expression memory is not<br />imposed. It emerges.</>,
  "chart_memory.png",
  <>Correlation between relatives collapses once the promoter switches faster than the cell divides.</>);

const Persist = chartScene("persist", "drug tolerance", ORANGE,
  <>One resistance gene reproduces<br />persister biology.</>,
  "chart_kill.png",
  <>Biphasic kill curves, and sisters share the fate of their lineage 4.8&times; more often than chance.</>);

const Calib = chartScene("calib", "calibration", BLUE,
  <>Fitted to published<br />time-lapse data.</>,
  "chart_fates.png",
  <>U2OS cells under cisplatin. Bars are the model, ticks are the measurement, one concentration held out.</>);

const Ident = chartScene("ident", "honest inference", PURPLE,
  <>And it reports what the<br />data cannot determine.</>,
  "chart_ident.png",
  <>Every outlined combination of memory and resistant fraction fits the same data within measurement noise.</>,
  470);

const Clinic = chartScene("clinic", "clinical schedules", ORANGE,
  <>Continuous or intermittent?<br />It depends on the cells.</>,
  "chart_schedules.png",
  <>Reproduces the melanoma trial outcome once you set how resistant cells grow with and without drug.</>);

const Close: React.FC = () => {
  const f = useCurrentFrame();
  return (
    <Scene durationInFrames={D.close}>
      <Reveal>
        <div style={{ fontSize: 76, fontWeight: 800, color: INK, letterSpacing: -1.8 }}>
          Biomodelling.jl <span style={{ color: ORANGE }}>2.0</span>
        </div>
      </Reveal>
      <Reveal delay={8}>
        <Sub size={34}>Mechanistic single-cell simulation, calibration and schedule optimisation &mdash; in Julia.</Sub>
      </Reveal>
      <Reveal delay={20} style={{ marginTop: 40 }}>
        <div style={{ display: "flex", gap: 14, flexWrap: "wrap" }}>
          {["open source", "MIT licence", "calibrated to published data", "white paper included"].map((c) => (
            <div key={c} style={{
              fontSize: 25, fontWeight: 600, color: INK2, border: "1.6px solid #e0ded8",
              borderRadius: 999, padding: "11px 22px", background: "#fff",
            }}>{c}</div>
          ))}
        </div>
      </Reveal>
      <Reveal delay={34} style={{ marginTop: 46 }}>
        <div style={{ fontSize: 31, fontWeight: 600, color: INK, opacity: Math.min(1, 0.55 + f / 90) }}>
          github.com/ayoublasri/Biomodelling.jl
        </div>
        <div style={{ fontSize: 26, fontWeight: 400, color: INK2, marginTop: 12 }}>
          Ayoub Lasri
        </div>
      </Reveal>
    </Scene>
  );
};

export const Showcase: React.FC = () => (
  <AbsoluteFill style={{ background: PAPER }}>
    <Fonts />
    <Sequence from={starts.hook} durationInFrames={D.hook}><Hook /></Sequence>
    <Sequence from={starts.gap} durationInFrames={D.gap}><Gap /></Sequence>
    <Sequence from={starts.mech} durationInFrames={D.mech}><Mech /></Sequence>
    <Sequence from={starts.memory} durationInFrames={D.memory}><Memory /></Sequence>
    <Sequence from={starts.persist} durationInFrames={D.persist}><Persist /></Sequence>
    <Sequence from={starts.calib} durationInFrames={D.calib}><Calib /></Sequence>
    <Sequence from={starts.ident} durationInFrames={D.ident}><Ident /></Sequence>
    <Sequence from={starts.clinic} durationInFrames={D.clinic}><Clinic /></Sequence>
    <Sequence from={starts.close} durationInFrames={D.close}><Close /></Sequence>
    <Progress total={TOTAL} />
  </AbsoluteFill>
);
