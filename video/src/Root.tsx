import React from "react";
import { Composition } from "remotion";
import { Showcase, TOTAL } from "./Showcase";

export const RemotionRoot: React.FC = () => (
  <>
    <Composition
      id="Showcase"
      component={Showcase}
      durationInFrames={TOTAL}
      fps={30}
      width={1080}
      height={1080}
    />
  </>
);
