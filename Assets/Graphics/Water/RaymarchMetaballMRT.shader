Shader "Custom/RaymarchMetaballMRT"
{
    Properties
    {
        _MainTex ("Texture", 2D) = "white" {}
        _Color ("Color", Color) = (1, 1, 1, 1)
        _Cube("Reflection Map", Cube) = ""{}
        _Radius("Radius", float) = 0.5
        _K("K", float) = 0.2
    }
    SubShader
    {
        Tags { "RenderType"="Transparent" "Queue"="Transparent" }
        Blend SrcAlpha OneMinusSrcAlpha
        ZTest Always
        ZWrite OFF
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma target 3.0

            #include "UnityCG.cginc"

            #define MAX_STEPS 100
            #define MAX_DIST 50
            #define SURF_DIST 0.001

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
                float3 ro : TEXCOORD1;
                float3 hitPos : TEXCOORD2;
                float3 light : TEXCOORD3;

                float3x3 viewMatrixIT : TEXCOORD4;
            };

            sampler2D _MainTex;
            samplerCUBE _Cube;
            float4 _MainTex_ST;
            sampler2D _CameraDepthTexture;
            float _Radius;
            float _K;

            StructuredBuffer<float3> _ParticlePoses;
            int _ParticleCount = 0;

            float4 _Color;

            float3x3 Inverse3x3(float3x3 m) {
                // 행렬식 계산
                float det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

                // 역행렬이 존재하는지 확인
                if (abs(det) < 1e-6) {
                    return float3x3(1, 0, 0, 0, 1, 0, 0, 0, 1); // 단위 행렬 반환 (역행렬이 없을 때)
                }

                // 여인수 전치 행렬 계산
                float3x3 adj;
                adj[0][0] =  (m[1][1] * m[2][2] - m[2][1] * m[1][2]);
                adj[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]);
                adj[0][2] =  (m[0][1] * m[1][2] - m[0][2] * m[1][1]);

                adj[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]);
                adj[1][1] =  (m[0][0] * m[2][2] - m[0][2] * m[2][0]);
                adj[1][2] = -(m[0][0] * m[1][2] - m[1][0] * m[0][2]);

                adj[2][0] =  (m[1][0] * m[2][1] - m[2][0] * m[1][1]);
                adj[2][1] = -(m[0][0] * m[2][1] - m[2][0] * m[0][1]);
                adj[2][2] =  (m[0][0] * m[1][1] - m[1][0] * m[0][1]);

                // 역행렬 계산: (1 / det) * adj
                return (1.0 / det) * adj;
            }


            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = TRANSFORM_TEX(v.uv, _MainTex);

                // calculate raymarch in world coord
                o.ro = float4(_WorldSpaceCameraPos, 1);
                o.hitPos = mul(unity_ObjectToWorld, v.vertex);
                o.light = _WorldSpaceLightPos0;

                // calculate view matrix's inverse transpos
                float3x3 viewMatrix = (float3x3) UNITY_MATRIX_V;
                o.viewMatrixIT = transpose(Inverse3x3(viewMatrix));
                return o;
            }

            float Metaball(float3 d, float3 c, float r) {
                return r / (length(d-c));   
            }

            // exponential
            float smin( float a, float b, float k )
            {
                k *= 1.0;
                float r = exp2(-a/k) + exp2(-b/k);
                return -k*log2(r);
            }

            float GetDist(float3 p) {
                float d = 0;
                {
                    float d1 = length(p - _ParticlePoses[0]) - _Radius;
                    float d2 = length(p - _ParticlePoses[1]) - _Radius;
                    d = smin(d1, d2, _K);
                }
                for (int i = 2; i < _ParticleCount; i++) {
                    float d1 = length(p - _ParticlePoses[i]) - _Radius;
                    d = smin(d, d1, _K);
                }
                return d;
            }

            float Raymarch(float3 ro, float3 rd) {
                float dO = 0;
                float dS;
                for (int i = 0; i < MAX_STEPS; i++) {
                    float3 p = ro + dO * rd;
                    dS = GetDist(p);
                    dO += dS;

                    if (dS < SURF_DIST || dO > MAX_DIST) break;
                }
                return dO;
            }

            float3 GetNormal(float3 p) {
                float2 e = float2(1e-2, 0);
                float3 n = GetDist(p) - float3(
                    GetDist(p - e.xyy),
                    GetDist(p - e.yxy),
                    GetDist(p - e.yyx)
                    );
                return normalize(n);
            }

            void frag (v2f i, out float4 color : SV_Target0, out float4 normal : SV_Target1)
            {
                float2 uv = i.uv - 0.5;
                float3 ro = i.ro;
                float3 rd = normalize(i.hitPos - ro);

                half d = Raymarch(ro, rd);
                fixed4 tex = tex2D(_MainTex, i.uv);
                fixed4 col = 0;
                half m = dot(uv, uv);

                if (d >= MAX_DIST) {
                    discard;
                } else {
                    half3 p = ro + rd * d;
                    float3 n = GetNormal(p);

                    // depth test
                    float4 view = mul(UNITY_MATRIX_V, float4(p, 1));
                    float4 proj = mul(UNITY_MATRIX_P, view);
                    proj /= proj.w;

                    half sceneDepth = tex2D(_CameraDepthTexture, float2(0.5 + proj.x / 2, 0.5 - proj.y / 2));
                    if (proj.z < sceneDepth)
                        discard;

                    // color
                    col.rgb = _Color;
                    col.a = 0.1;
                    
                    // fresnel
                    half cosNR = saturate(dot(n, normalize(-rd)));
                    half r0 = pow((1.0003 - 1.33) / (1.0003 + 1.33), 2);
                    half r = r0 + (1 - r0) * pow(1 - cosNR, 2);
                    
                    // environment
                    float3 refl = reflect(rd, n);
                    half4 ecol = UNITY_SAMPLE_TEXCUBE(unity_SpecCube0, refl); 
                    half3 finalEColor = DecodeHDR(ecol, unity_SpecCube0_HDR);
                    //finalEColor = texCUBE(_Cube, n);

                    col.a += r;
                    col.rgb = finalEColor * r;

                    // specular
                    half3 h = -normalize(rd - i.light);
                    half sp = pow(saturate(dot(h, n)), 1000);
                    col += sp;
                    
                    color = col;
                    half3 n_cam = mul(i.viewMatrixIT, n);
                    n_cam.xyz = normalize(n_cam.xyz);
                    n_cam.xyz += 1.0f;
                    n_cam.xyz /= 2.0f;
                    
                    normal = float4(n_cam, 1);
                    //color = half4(n_cam.xyz, 1);
                } 
            }

            
            // https://www.youtube.com/watch?v=S8AWd66hoCo&ab_channel=TheArtofCode
            ENDCG
        }
    }
}
