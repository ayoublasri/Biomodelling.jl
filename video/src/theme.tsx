import React from "react";
import { interpolate, spring, staticFile, useCurrentFrame, useVideoConfig } from "remotion";

export const PAPER = "#fcfbf9";
export const INK = "#0b0b0b";
export const INK2 = "#575652";
export const LINE = "#e4e2dc";
export const BLUE = "#2a78d6";
export const ORANGE = "#eb6834";
export const GREEN = "#1baf7a";
export const RED = "#e34948";
export const PURPLE = "#4a3aa7";

export const Fonts: React.FC = () => (
  <style>{`
    @font-face { font-family: 'Inter'; font-weight: 400;
      src: url('${staticFile("fonts/Inter-400.ttf")}') format('truetype'); }
    @font-face { font-family: 'Inter'; font-weight: 600;
      src: url('${staticFile("fonts/Inter-600.ttf")}') format('truetype'); }
    @font-face { font-family: 'Inter'; font-weight: 800;
      src: url('${staticFile("fonts/Inter-800.ttf")}') format('truetype'); }
  `}</style>
);

/** Fade + rise, held after it lands. `delay` is in frames. */
export const Reveal: React.FC<{
  delay?: number; children: React.ReactNode; y?: number; style?: React.CSSProperties;
}> = ({ delay = 0, children, y = 26, style }) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();
  const s = spring({ frame: frame - delay, fps, config: { damping: 200, mass: 0.6 }, durationInFrames: 26 });
  return (
    <div style={{ opacity: s, transform: `translateY(${(1 - s) * y}px)`, ...style }}>{children}</div>
  );
};

/** Wraps a scene: paper background, padding, and a crossfade at both ends. */
export const Scene: React.FC<{ durationInFrames: number; children: React.ReactNode }> = ({
  durationInFrames, children,
}) => {
  const frame = useCurrentFrame();
  const opacity = interpolate(
    frame,
    [0, 8, durationInFrames - 9, durationInFrames - 1],
    [0, 1, 1, 0],
    { extrapolateLeft: "clamp", extrapolateRight: "clamp" }
  );
  return (
    <div style={{
      position: "absolute", inset: 0, opacity, background: PAPER,
      fontFamily: "Inter, sans-serif", display: "flex", flexDirection: "column",
      justifyContent: "center", padding: "0 86px", boxSizing: "border-box",
    }}>{children}</div>
  );
};

export const Kicker: React.FC<{ children: React.ReactNode; color?: string }> = ({ children, color = ORANGE }) => (
  <div style={{
    fontSize: 25, fontWeight: 600, letterSpacing: 2.6, textTransform: "uppercase",
    color, marginBottom: 22,
  }}>{children}</div>
);

export const Headline: React.FC<{ children: React.ReactNode; size?: number }> = ({ children, size = 66 }) => (
  <div style={{ fontSize: size, fontWeight: 800, lineHeight: 1.13, color: INK, letterSpacing: -1.4 }}>
    {children}
  </div>
);

export const Sub: React.FC<{ children: React.ReactNode; size?: number }> = ({ children, size = 33 }) => (
  <div style={{ fontSize: size, fontWeight: 400, lineHeight: 1.42, color: INK2, marginTop: 26 }}>
    {children}
  </div>
);

/** A chart image on a soft card. */
export const ChartCard: React.FC<{ src: string; delay?: number; maxHeight?: number }> = ({
  src, delay = 0, maxHeight = 500,
}) => {
  const frame = useCurrentFrame();
  const { fps } = useVideoConfig();
  const s = spring({ frame: frame - delay, fps, config: { damping: 200, mass: 0.7 }, durationInFrames: 30 });
  return (
    <div style={{
      opacity: s, transform: `scale(${0.965 + s * 0.035})`, marginTop: 34,
      display: "flex", justifyContent: "center",
    }}>
      <img src={staticFile(src)} style={{
        maxWidth: "100%", maxHeight, objectFit: "contain",
        borderRadius: 16, border: `1px solid ${LINE}`, background: "#fff",
        boxShadow: "0 12px 34px rgba(18,16,12,0.07)",
      }} />
    </div>
  );
};

/** Thin progress bar across the bottom of every frame. */
export const Progress: React.FC<{ total: number }> = ({ total }) => {
  const frame = useCurrentFrame();
  const p = Math.min(frame / total, 1);
  return (
    <div style={{ position: "absolute", left: 0, right: 0, bottom: 0, height: 5, background: "#eeece7" }}>
      <div style={{ width: `${p * 100}%`, height: "100%", background: ORANGE }} />
    </div>
  );
};
