import React from "react";
import { interpolate, useCurrentFrame } from "remotion";
import { BLUE, GREEN, INK2, LINE, ORANGE, RED } from "./theme";

/** Stable pseudo-random molecule positions (no RNG: the same frame must always draw the same picture). */
const dots = (n: number, seed: number) =>
  Array.from({ length: n }, (_, i) => {
    const a = Math.sin((i + 1) * 12.9898 + seed) * 43758.5453;
    const b = Math.sin((i + 1) * 78.233 + seed) * 12345.6789;
    const r = Math.sqrt(Math.abs(a % 1)) * 0.82;
    const th = (Math.abs(b % 1)) * Math.PI * 2;
    return { x: Math.cos(th) * r, y: Math.sin(th) * r };
  });

/** Scene 1: two identical cells, one drug, one survivor. */
export const TwoCells: React.FC<{ duration: number }> = ({ duration }) => {
  const f = useCurrentFrame();
  const t = f / duration;
  const appear = interpolate(t, [0, 0.16], [0, 1], { extrapolateRight: "clamp" });
  const wash = interpolate(t, [0.34, 0.5], [0, 1], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
  const kill = interpolate(t, [0.54, 0.74], [0, 1], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
  const pulse = 1 + 0.035 * Math.sin(t * Math.PI * 7);
  const R = 96;
  const cell = (cx: number, dying: boolean, protectedCell: boolean) => (
    <g opacity={appear} transform={`translate(${cx},170)`}>
      <circle r={R} fill="#fff" stroke={dying ? RED : INK2}
        strokeWidth={dying ? 3 : 2.4} opacity={dying ? 1 - kill * 0.75 : 1}
        transform={`scale(${dying ? 1 - kill * 0.22 : protectedCell ? pulse : 1})`} />
      {protectedCell && (
        <circle r={R + 13} fill="none" stroke={GREEN} strokeWidth={3.4}
          opacity={kill * 0.85} strokeDasharray="7 9" />
      )}
      {dots(16, dying ? 3 : 9).map((d, i) => (
        <circle key={i} cx={d.x * R} cy={d.y * R} r={6.4}
          fill={i < 4 ? BLUE : "#b9b7b1"} opacity={dying ? 1 - kill : 1} />
      ))}
      {dying && (
        <g opacity={kill} stroke={RED} strokeWidth={8} strokeLinecap="round">
          <line x1={-34} y1={-34} x2={34} y2={34} />
          <line x1={34} y1={-34} x2={-34} y2={34} />
        </g>
      )}
    </g>
  );
  return (
    <svg viewBox="0 0 900 340" style={{ width: "100%", height: 340 }}>
      <rect x={0} y={0} width={900} height={340 * wash} fill={ORANGE} opacity={0.085} />
      {wash > 0.02 && (
        <text x={872} y={38} textAnchor="end" fontSize={25} fontWeight={600} fill={ORANGE} opacity={wash}>
          drug
        </text>
      )}
      {cell(272, true, false)}
      {cell(628, false, true)}
    </svg>
  );
};

/** Scene 3: grow, replicate, divide, and both daughters inherit the promoter state. */
export const CellCycle: React.FC<{ duration: number }> = ({ duration }) => {
  const f = useCurrentFrame();
  const t = f / duration;
  const grow = interpolate(t, [0.06, 0.40], [1, 1.34], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
  const rep = interpolate(t, [0.40, 0.55], [0, 1], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
  const split = interpolate(t, [0.60, 0.86], [0, 1], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
  const nmol = Math.round(interpolate(t, [0.06, 0.58], [9, 22], { extrapolateLeft: "clamp", extrapolateRight: "clamp" }));
  const R = 78;
  const sep = split * 205;
  // volume halves at division, so the radius falls by 2^(-1/3); interpolated so the split does not pop
  const dScale = grow * (1 - 0.206 * split);
  const labels: [string, number][] = [["grow", 0.12], ["replicate", 0.42], ["divide", 0.63], ["inherit", 0.86]];
  const daughter = (dx: number, seed: number) => (
    <g transform={`translate(${450 + dx},168) scale(${dScale})`}>
      <circle r={R} fill="#fff" stroke={INK2} strokeWidth={2.4} />
      <circle cx={0} cy={-R * 0.52} r={11} fill={BLUE} />
      {dots(Math.max(4, Math.round(nmol / (1 + split))), seed).map((d, i) => (
        <circle key={i} cx={d.x * R} cy={d.y * R * 0.92} r={6.2} fill={i % 3 === 0 ? ORANGE : "#b9b7b1"} />
      ))}
    </g>
  );
  return (
    <svg viewBox="0 0 900 340" style={{ width: "100%", height: 340 }}>
      <line x1={70} y1={296} x2={830} y2={296} stroke={LINE} strokeWidth={3} />
      {labels.map(([txt, at], i) => {
        const on = interpolate(t, [at, at + 0.07], [0, 1], { extrapolateLeft: "clamp", extrapolateRight: "clamp" });
        return (
          <g key={i} opacity={on}>
            <circle cx={128 + i * 215} cy={296} r={8} fill={ORANGE} />
            <text x={128 + i * 215} y={332} textAnchor="middle" fontSize={25} fontWeight={600} fill={INK2}>
              {txt}
            </text>
          </g>
        );
      })}
      {split > 0.02 ? (
        <>{daughter(-sep, 9)}{daughter(sep, 4)}</>
      ) : (
        <g transform={`translate(450,168) scale(${grow})`}>
          <circle r={R} fill="#fff" stroke={INK2} strokeWidth={2.4} />
          <circle cx={rep > 0.5 ? -13 : 0} cy={-R * 0.52} r={11} fill={BLUE} />
          {rep > 0.5 && <circle cx={13} cy={-R * 0.52} r={11} fill={BLUE} opacity={rep} />}
          {dots(nmol, 9).map((d, i) => (
            <circle key={i} cx={d.x * R} cy={d.y * R * 0.92} r={6.2} fill={i % 3 === 0 ? ORANGE : "#b9b7b1"} />
          ))}
        </g>
      )}
    </svg>
  );
};
